#include "pasm.h"

AlgoPASM::AlgoPASM(EntityId modelID,size_t nface, const int* faces)
{
  _lastError = MR_Enums::eSucceed;
  _maxId = 0;

  auto modePtr = GetModel(modelID);
  _modePtr = &(*modePtr);
  if (!_modePtr) {
    _lastError = MR_Enums::eInvalidInput;
    return;
  }

  _faceMesh.reserve(nface);
  IntSet faceSet;
  IntSet edgeSet;
  for (size_t i = 0; i < nface; i++) {
    if (faceSet.insert(faces[i]).second) {
      TpFace face;
      if (_modePtr->getEntity(faces[i], face)) {
        MeshDataPtr data(new MeshData);
        data->face = face;
        _faceMesh.push_back(data);
        TpFaceEdgeIter faceEdge(face);
        TpEdge edge;

        for (faceEdge.first(); faceEdge.next(edge);) {
          if (edgeSet.insert(edge.getIndex()).second) {
            _edgeMesh[edge.getIndex()].edge = edge;
            TpVertex vx[3];
            edge.getVertex(vx);
            for (int i = 0; i < 2; i++) {
              if (vx[i].isValid()) {
                _maxId = std::max(_maxId, vx[i].getIndex());
              }
            }
          }
        }
      }
    }
  }

  for (auto it : _faceMesh) {
    _faceMeshMap[it->face.getIndex()] = it;
  }

}

void AlgoPASM::setup(const char* ctrlstr)
{
  double pmin[3], pmax[3];
  _modePtr->bbox(pmin, pmax);
  BBox3D box;
  box += pmin;
  box += pmax;
  _gamma = box.size(true);

  _dfMax = _gamma * 0.1;
  _dfMin = _dfMax * 0.001;
  _angle = 20.0;
  _angleDet = cos(ToRadian(_angle* 0.25));
  _alpha_c = 0.2;
  _alpha_nn = 0.1;
  _alpha_ne = 0.05;

  _mergedNodeId = -1;

  _d_c_relative = 0.2;
  _d_c = 0.0;
  _d_m = 0.0;
  _alpha_c = 0.2;

  _dfMaxBKG = _gamma * 0.05;
  _eps = _dfMaxBKG * 0.001;
  _dfMinBKG = _dfMaxBKG * 0.001;

  Json::Value ctrl;
  if (ctrlstr) {
    Json::Reader reader;
    reader.parse(ctrlstr, ctrl);
  }
  else {
    return;
  }

  if (ctrl.empty())
    return;

    double tmp;
    if (tGetValue(ctrl, "meshSizeMax", tmp) && tmp > Gloabals::ZERO) {
      _dfMax = tmp;
    }
    else if (tGetValue(ctrl, "meshSizePercentOfBBox", tmp) && tmp > Gloabals::ZERO && tmp < 99.0) {
      _dfMax = _gamma*tmp / 100.0;
    }

    if (tGetValue(ctrl, "meshSizeMin", tmp) && tmp > Gloabals::ZERO) {
      _dfMin = tmp;
    }
    else if (tGetValue(ctrl, "meshSizeMinPercentOfMax", tmp) && tmp > Gloabals::ZERO && tmp < 99.0) {
      _dfMin = _dfMax * tmp / 100.0;
    }
    else {
      _dfMin = _dfMax * 0.001;
    }

    if (tGetValue(ctrl, "turnAngle", tmp) && tmp > Gloabals::ZERO) {
      double angle = tmp;
      if (angle <= 0 || angle >= 90)
        angle = 20;
      _angle = angle;
      _angleDet = cos(ToRadian(_angle*0.25));
    }

    if (tGetValue(ctrl, "alpha_c", tmp) && tmp > Gloabals::ZERO) {
      _alpha_c = tmp;
    }

    if (tGetValue(ctrl, "alpha_nn", tmp) && tmp > Gloabals::ZERO) {
      _alpha_nn = tmp;
    }

    if (tGetValue(ctrl, "alpha_ne", tmp) && tmp > Gloabals::ZERO) {
      _alpha_ne = tmp;
    }

    if (tGetValue(ctrl, "d_c", tmp) && tmp > Gloabals::ZERO) {
      _d_c = tmp;
      _d_c_relative = 0.0;
    }
    else if (tGetValue(ctrl, "d_c_relative", tmp) && tmp > Gloabals::ZERO) {
      _d_c = 0.0;
      _d_c_relative = tmp;
    }

    if (tGetValue(ctrl, "d_s", tmp) && tmp > Gloabals::ZERO) {
      _d_s = tmp;
    }

    if (tGetValue(ctrl, "d_m", tmp) && tmp > Gloabals::ZERO) {
      _d_m = tmp;
    }
    else {
      _d_m = _d_s;
    }

}

int AlgoPASM::operator()(const char* ctrlstr)
{
  setup(ctrlstr);
  aligned_curve_meshing();
  aligned_surface_meshing();

  return _lastError;
}

void AlgoPASM::setup_acm()
{
  auto lmdCollectBKGCurveMesh = [this]() {
    //collect vertex nodes.
    TpVertex vx[2];
    int maxVertId = -1;
    Real5 anode;
    for (auto& it : _curvedataarray) {
      it.edge.getVertex(vx);
      for (const auto& vt : vx) {
        int vid = vt.getIndex();
        maxVertId = std::max(maxVertId, vid);
        _validVertexNodeSet.insert(vid);
        if (vid + 1 > _curvenodes.size())
          _curvenodes.resize(vid + 1);
        vt.getPosition(anode.data());
        _curvenodes[vid] = anode;
      }
    }

    // collect curve nodes and segments.
    const int ncurve = _curvedataarray.size();
    for (int ci = 0; ci < ncurve; ++ci) {
      CurveData& it = _curvedataarray[ci];
      it.edge.getVertex(vx);
      if (!it.edge.isForward()) std::swap(vx[0], vx[1]);
      int vids[2] = { vx[0].getIndex(),vx[1].getIndex() };
      assert(it.tarray0.size() > 1);
      DIntArray nodeids;
      nodeids.push_back(vids[0]);
      int nnode = it.tarray0.size();
      for (int i = 1; i < nnode - 1; ++i) {
        it.edge.eval(it.tarray0[i], anode.data());
        _curvenodes.push_back(anode);
        nodeids.push_back(_curvenodes.size() - 1);
      }
      nodeids.push_back(vids[1]);

      int nid = nodeids.size();
      for (int i = 0; i < nid - 1; ++i) {
        _segments.push_back(Int3(nodeids[i], nodeids[i + 1], ci));
      }

    }

  };

  auto lmdSetupBKGCurveMeshSize = [this]() {

    DIntArray curvenodecnt(_curvenodes.size(), 0);
    for (auto& nd : _curvenodes) {
      nd[3] = 0.0;
      nd[4] = Gloabals::MAX_REAL;
    }

    for (const auto& seg : _segments) {
      Point3D pt0(_curvenodes[seg[0]].data());
      Point3D pt1(_curvenodes[seg[1]].data());
      double len = pt0.distance(pt1);
      if (len < Gloabals::ZERO) {
        Messager::Msg(Messager::eWarning, "Degenerated curve %d found!", _curvedataarray[seg[2]].edge.getIndex());
        len = _eps;
      }
      for (int i = 0; i < 2; ++i) {
        _curvenodes[seg[i]][3] += len;
        ++curvenodecnt[seg[i]];
        _curvenodes[seg[i]][4] = std::min(_curvenodes[seg[i]][4], len);
      }
    }

    int nnode = _curvenodes.size();
    for (int i = 0; i < nnode; ++i) {
      if (curvenodecnt[i] < 1) 
        continue;
      _curvenodes[i][3] /= (double)curvenodecnt[i];
    }

  };

  // do initial curve division.
  for (auto& it : _edgeMesh) {
    _curvedataarray.push_back(CurveData(it.second.edge));
  }

  const int ncurve = _curvedataarray.size();
#pragma omp parallel for
  for (int ci = 0; ci < ncurve; ++ci) {
    auto& curvedata = _curvedataarray[ci];
    TpEdgeWrapper myedge(curvedata.edge);
    double intv[2];
    curvedata.edge.interval(intv);
    CurveDivide<TpEdgeWrapper, Point3D>(myedge, _dfMaxBKG, _dfMinBKG, cos(ToRadian(_angle4background* 0.25)),
      intv, curvedata.tarray4Itsct, false);
    _curveMeshFun(myedge, intv, curvedata.tarray0);
  }

  lmdCollectBKGCurveMesh();
  lmdSetupBKGCurveMeshSize();

  for (auto& cdata : _curvedataarray) {
    MRAssert(cdata.tarray1.size() == 2);
    cdata.tarray1[0]._df = _curvenodes[cdata.tarray1[0]._id][4];
    cdata.tarray1[1]._df = _curvenodes[cdata.tarray1[1]._id][4];
  }

  _CNodeMaxId = _maxId;
}

void AlgoPASM::aligned_curve_meshing()
{
  setup_acm();

  tBoxHash<double, 3> segmenthash(_segments.size(),
    [this](size_t i, tBoxHash<double, 3>::Box& box) {
    auto seg = _segments[i];
    box.reset();
    box |= _curvenodes[seg[0]];
    box |= _curvenodes[seg[1]];
    return true;
  });

  synchronize_curves_for_vertices(segmenthash);
  insert_touching_points(segmenthash);
  synchronize_curves_for_inner_nodes(segmenthash);
  confirm_curve_mesh();
}

void AlgoPASM::confirm_curve_mesh() 
{
  auto lmdConfirmCurveMesh = [this](CurveData& curvedata)
  {
    TpVertex vx[2];
    EdgeData edgeData;
    curvedata.edge.getVertex(vx);
    if (!curvedata.edge.isForward()) {
      std::swap(vx[0], vx[1]);
    }
    edgeData.edge = curvedata.edge;
    auto& pts = curvedata.tarray1;
    edgeData.pts.resize(pts.size());
    edgeData.pts[0].first = vx[0].getIndex();
    edgeData.pts[0].second = pts[0]._t;
    auto sz = pts.size() - 1;
    for (size_t j = 1; j < sz; j++) {
      edgeData.pts[j].first = ++_maxId;
      edgeData.pts[j].second = pts[j]._t;
    }
    edgeData.pts[sz].first = vx[1].getIndex();
    edgeData.pts[sz].second = pts.back()._t;
    curvedata.edge.check(TpEntity::eMeshDone, true);
    _edgeMesh[curvedata.edge.getIndex()] = edgeData;
  };

  int maxIdBKG = _maxId;
  auto lmdConfirmBKGCurveMesh = [this, &maxIdBKG](const CurveData& curvedata)
  {
    TpVertex vx[2];
    EdgeData edgeData;
    auto& pts = curvedata.tarray4Itsct;
    edgeData.pts.resize(pts.size());
    edgeData.pts[0].first = vx[0].getIndex();
    edgeData.pts[0].second = pts[0];
    auto sz = pts.size() - 1;
    for (size_t j = 1; j < sz; j++) {
      edgeData.pts[j].first = ++maxIdBKG;
      edgeData.pts[j].second = pts[j];
    }
    edgeData.pts[sz].first = vx[1].getIndex();
    edgeData.pts[sz].second = pts.back();
    _edgeBKGMesh[curvedata.edge.getIndex()] = edgeData;
  };

  for (auto& curvedata : _curvedataarray) {
    lmdConfirmCurveMesh(curvedata);
    lmdConfirmBKGCurveMesh(curvedata);
  }

}

void AlgoPASM::synchronize_curves_for_inner_nodes(const tBoxHash<double, 3>& segmenthash)
{
  auto lmdDoCurveMesh = [this](CNodeArray& msh,
    TpEdgeWrapper& myedge,
    CNodeArray& nodesnew) {

    auto lmdSetCNodeArrayMeshSize = [&myedge](CNodeArray& msh) {
      int nnode = msh.size();
      double len;
      Point3D pt0, pt1;
      for (int ni = 0; ni < nnode - 1; ++ni) {
        pt0 = myedge.eval(msh[ni]._t);
        pt1 = myedge.eval(msh[ni + 1]._t);
        len = pt0.distance(pt1);
        msh[ni]._df = std::min(msh[ni]._df, len);
        msh[ni + 1]._df = std::min(msh[ni + 1]._df, len);
        MRAssert(msh[ni]._df > Gloabals::ZERO);
        MRAssert(msh[ni + 1]._df > Gloabals::ZERO);
      }
    };

    DIntArray nodesnewid;
    CNodeArray msh1;
    DRealArray tmp;
    CNode cnode;
    cnode._df = Gloabals::MAX_REAL;
    int sz = msh.size();
    for (int i = 0; i < sz - 1; ++i) {
      msh1.push_back(msh[i]);
      tmp.clear();
      _curveMeshFun(myedge, Real2(msh[i]._t, msh[i + 1]._t), tmp);
      int sztmp = tmp.size();
      for (int ti = 1; ti < sztmp - 1; ++ti) {
        cnode._id = ++_CNodeMaxId;
        cnode._t = tmp[ti];
        msh1.push_back(cnode);
        nodesnewid.push_back(msh1.size() - 1);
      }
    }
    cnode._id = ++_CNodeMaxId;
    cnode._t = tmp.back();
    msh1.push_back(cnode);
    lmdSetCNodeArrayMeshSize(msh1);
    msh.swap(msh1);

    for (int ni : nodesnewid) {
      nodesnew.push_back(msh[ni]);
    }
  };

  CNodeArray nodesnew;
  auto lmdDoMesh = [&](int ci) {
    CurveData& curvedata = _curvedataarray[ci];
    if (curvedata.isMeshed) 
      return;
    TpEdgeWrapper myedge(curvedata.edge);
    nodesnew.clear();
    lmdDoCurveMesh(curvedata.tarray1, myedge, nodesnew);
    do_synchronize_curves_around_nodes(segmenthash, nodesnew, ci);
    curvedata.isMeshed = true;
  };

  // sort.
  const int ncurvedata = _curvedataarray.size();
  std::vector<std::pair<int, double>> accessor;
  accessor.reserve(ncurvedata);
  for (int i = 0; i < ncurvedata; ++i) {
    accessor.push_back(std::make_pair(i, _curvedataarray[i].edge.length()));
  }

  std::sort(accessor.begin(), accessor.end(),
    [](const std::pair<int, double>& a, const std::pair<int, double>& b) {return a.second < b.second; });

  for (int ci = 0; ci < ncurvedata; ++ci) {
    lmdDoMesh(accessor[ci].first);
  }// for (int ci = 0; ci < ncurvedata; ++ci)

}

void AlgoPASM::synchronize_curves_for_vertices(const tBoxHash<double, 3>& segmenthash)
{
  std::vector<NodeInfo> nodesinfo;
  nodesinfo.reserve(_validVertexNodeSet.size());
  NodeInfo ninfo;
  for (auto& myvertidx : _validVertexNodeSet) {
    Real5& coord = _curvenodes[myvertidx];
    IntSet ciset;
    CNode cnd;
    cnd._df = coord[3];
    MRAssert(cnd._df > Gloabals::ZERO);
    cnd._id = myvertidx;
    ninfo.cnd = cnd;
    ninfo.cset = std::make_shared<IntSet>();
    ninfo.pt3d = Point3D(coord.data());
    ninfo.tangent = Vector3D(.0, .0, .0);
    ninfo.filter = [this, myvertidx](int i)->bool {
      TpEdge& edge = _curvedataarray[_segments[i][2]].edge;
      TpVertex vx[2];
      edge.getVertex(vx);
      return (vx[0].getIndex() != myvertidx && vx[1].getIndex() != myvertidx);
    };
    nodesinfo.push_back(ninfo);
  }
  synchronize_curves_around_nodes(segmenthash, nodesinfo);

}

void AlgoPASM::synchronize_curves_around_nodes(const tBoxHash<double, 3>& segmenthash, std::vector<NodeInfo>& nodesinfo)
{
  if (nodesinfo.empty())
    return;

  const int maxIter = 3;
  int cnt = 0;
  std::vector<NodeInfo> nodesinfo2;
  for (;;) {
    int ninfo = nodesinfo.size();
#pragma omp parallel for
    for (int ni = 0; ni < ninfo; ++ni) {
      synchronize_curves_around_node(segmenthash, nodesinfo[ni].pt3d,
        nodesinfo[ni].tangent,
        nodesinfo[ni].cnd,
        nodesinfo[ni].cset,
        nodesinfo[ni].filter,
        nodesinfo2);
    }
    if (nodesinfo2.empty())
      break;
    if (cnt++ > maxIter)
      break;
    nodesinfo.swap(nodesinfo2);
    nodesinfo2.clear();
  }

}

void AlgoPASM::do_synchronize_curves_around_nodes(const tBoxHash<double, 3>& segmenthash, const CNodeArray& nodesnew, int ci)
{
  if (nodesnew.empty())
    return;

  CurveData& curvedata = _curvedataarray[ci];
  int nnode = nodesnew.size();
  std::vector<NodeInfo> nodesinfo;
  nodesinfo.reserve(nnode);
  NodeInfo ninfo;
  for (int ni = 0; ni < nnode; ++ni) {
    CNode cnd = nodesnew[ni];
    Point3D pt;
    curvedata.edge.eval(cnd._t, pt.data());
    Vector3D dir(.0, .0, .0);
    curvedata.edge.tangent(cnd._t, dir.data());
    ninfo.cnd = cnd;
    ninfo.cset = std::make_shared<IntSet>();
    ninfo.cset->insert(ci);
    ninfo.pt3d = pt;
    ninfo.tangent = dir;
    ninfo.filter = [ci, this](int i) {return _segments[i][2] != ci; };
    nodesinfo.push_back(ninfo);
  }
  synchronize_curves_around_nodes(segmenthash, nodesinfo);

}

double AlgoPASM::get_curve_mesh_size(CurveData& curvedata, double t)
{
  TpEdge& edge = curvedata.edge;
  DRealArray& tarray = curvedata.tarray0;
  auto pos = std::lower_bound(tarray.begin(), tarray.end(), t);
  Point3D pt0, pt1;
  if (pos == tarray.begin()) {
    edge.eval(*(pos + 1), pt0.data());
    edge.eval(*(pos), pt1.data());
  }
  else if (pos == tarray.end()) {
    edge.eval(*(pos - 2), pt0.data());
    edge.eval(*(pos - 1), pt1.data());
  }
  else {
    edge.eval(*(pos - 1), pt0.data());
    edge.eval(*(pos), pt1.data());
  }
  return pt0.distance(pt1);
};

bool AlgoPASM::insert_curve_node(CNode& tonedge, TpEdge& edge, CNodeArray& tarray, double rate, bool accurate)
{
  MRAssert(tonedge._df > 0.0);
  auto pos = std::lower_bound(tarray.begin(), tarray.end(), tonedge);
  if (pos == tarray.begin())
    return false;
  if (pos == tarray.end())
    return false;
  CNodeArray::iterator poses[2] = { pos,pos - 1 };
  Point3D ptonedge, pt1;
  edge.eval(tonedge._t, ptonedge.data());

  double dis2[2];
  for (int k = 0; k < 2; ++k) {
    if (accurate) {
      dis2[k] = edge.length(std::min(poses[k]->_t, tonedge._t), std::max(poses[k]->_t, tonedge._t));
    }
    else {
      edge.eval(poses[k]->_t, pt1.data());
      dis2[k] = ptonedge.distance2(pt1);
    }
    MRAssert(poses[k]->_df > Gloabals::ZERO);// valid mesh size ought to be provided.
    if (dis2[k] < poses[k]->_df*poses[k]->_df*rate*rate)
      return false;
  }
  tarray.insert(pos, tonedge);
  return true;
};

void AlgoPASM::synchronize_curves_around_node(const tBoxHash<double, 3>& segmenthash,
  const Point3D& pt,
  const Vector3D& dir,
  const CNode& cnd,
  std::shared_ptr<IntSet> ciset,
  std::function<bool(int)> filter,
  std::vector<NodeInfo>& nodesinfo2)
{

  auto lmdHasCloseSegment = [this, &pt](DIntArray& elems, double eps)->bool {
    // skip projection for distant curves.
    Point3D closestpt;
    double t;
    for (auto si : elems) {
      Real5 nd0 = _curvenodes[_segments[si][0]];
      Real5 nd1 = _curvenodes[_segments[si][1]];
      Segment3D seg(nd0.data(), nd1.data());
      Point3D closestpt = seg.closestPointTo(pt, t);
      if (t<-0.1 || t>1.1) 
        continue;
      if (closestpt.distance2(pt) > eps*eps) {
        continue;
      }
      else {
        return true;
      }
    }
    return false;
  };

  if (cnd._df < Gloabals::ZERO)
    return;
  double rad = get_curve_align_radius(cnd._df);
  if (rad < Gloabals::ZERO)
    return;
  tBoxHash<double, 3>::Box box;
  box |= pt.data();
  box.thicknessIt(rad);
  DIntArray includes;
  segmenthash.getBoxIncludes(box.min().data(), box.max().data(),
    [&includes, &filter](int i) {
    if (filter(i)) {
      includes.push_back(i);
    }
    return true;
  });
  if (includes.empty())
    return;

  bool checkdir = dir[0] + dir[0] + dir[0] > Gloabals::ZERO;
  // prepare data.
  boost::unordered_map<int, DIntArray> curvesegsmap;
  for (int si : includes) {
    int curveid = _segments[si][2];
    auto& elems = curvesegsmap[curveid];
    elems.push_back(si);
  }

  double rad2 = rad * rad;
  for (auto& it : curvesegsmap) {
    int dstci = it.first;
    CurveData& curvedata = _curvedataarray[dstci];

    if (!lmdHasCloseSegment(it.second, rad*3.0))
      continue;

    // do projection.
    Point3D ptonedge;
    double tonedge, dis;
    dis = curvedata.edge.project(pt.data(), ptonedge.data(), tonedge);
    if (dis > rad)
      continue;
    // check direction.
    if (checkdir && dis > _eps) {
      Vector3D vec(pt, ptonedge);
      vec.normalize();
      double dot = vec.dot(dir);
      if (dot > 0.414 || dot < -0.414)
        continue;
    }

    CNode nd;
    nd._df = get_curve_mesh_size(curvedata, tonedge);
    nd._t = tonedge;
    bool rt = false;
    {
      std::lock_guard<std::mutex> lock(_mutex);
      MRAssert(nd._df > Gloabals::ZERO);
      rt = insert_curve_node(nd, curvedata.edge, curvedata.tarray1, _alpha_c, false);
    }

    if (rt) {
      Vector3D dir1(.0, .0, .0);
      double intv[2];
      curvedata.edge.interval(intv);
      double eps = 0.001*std::abs(intv[1] - intv[0]);
      if (std::abs(tonedge - intv[0]) > eps&&std::abs(tonedge - intv[1]) > eps) {
        if (curvedata.edge.tangent(tonedge, dir1.data()) != MR_Enums::eSucceed)
          continue;
      }
      NodeInfo ninfo;
      ninfo.cnd = nd;
      ninfo.pt3d = ptonedge;
      ninfo.tangent = dir;
      ninfo.cset = ciset;
      ninfo.filter = [this, dstci](int i) {return _segments[i][2] != dstci; };
      std::lock_guard<std::mutex> lock(_mutex);
      nodesinfo2.push_back(ninfo);
    }

  }

};

void AlgoPASM::insert_touching_points(const tBoxHash<double, 3>& segmenthash)
{
  int ncurvedata = _curvedataarray.size();
#pragma omp parallel for
  for (int ci = 0; ci < ncurvedata; ++ci) {
    CurveData& curvedata = _curvedataarray[ci];
    TpEdgeWrapper myedge(curvedata.edge);
    int curveid = curvedata.edge.getIndex();
    DRealArray intersections;
    get_touching_points(segmenthash, myedge, ci, curvedata.tarray4Itsct, intersections);
    CNodeArray nodesnew;
    nodesnew.reserve(intersections.size());
    for (auto t : intersections) {
      CNode cnd;
      cnd._t = t;
      cnd._df = get_curve_mesh_size(curvedata, t);
      if (insert_curve_node(cnd, curvedata.edge, curvedata.tarray1, _alpha_c, false))
        nodesnew.push_back(cnd);
    }
    do_synchronize_curves_around_nodes(segmenthash, nodesnew, ci);
  }

}

struct ItsctInfo
{
  ItsctInfo() :dis(Gloabals::MAX_REAL), sc(-1.0) {}
  ItsctInfo(double a, double b) :dis(a), sc(b) {}

  bool operator < (const ItsctInfo& src) const
  {
    return this->dis < src.dis;
  }

  double dis;
  double sc;
};

void AlgoPASM::get_touching_points(const tBoxHash<double, 3>& segmenthash,
  TpEdgeWrapper& edge,
  int myci,
  const DRealArray& mytarray,
  DRealArray& intersections)
{
  boost::unordered_map<int, std::vector<ItsctInfo>> mindismap;

  auto lmdCollectMinDis = [&]() {

    double rad2 = _d_s * _d_s;
    double tol = _d_s * 1.5;
    double tol2 = tol * tol;
    double angtol = cos(ToRadian(10.0));
    Point3D pt0, pt1;
    tBoxHash<double, 3>::Box box;
    int nt = mytarray.size();
    for (int ti = 0; ti < nt - 1; ++ti) {
      Real2 intv(mytarray[ti], mytarray[ti + 1]);
      pt0 = edge.eval(intv[0]);
      pt1 = edge.eval(intv[1]);
      double len0 = pt0.distance2(pt1);
      if (len0 < Gloabals::ZERO) 
        continue;
      Segment3D seg0(pt0, pt1);
      Vector3D dir0(pt0, pt1);
      dir0 /= sqrt(len0);

      box.reset();
      box |= pt0;
      box |= pt1;
      box.thicknessIt(tol);
      DIntArray includes;
      segmenthash.getBoxIncludes(box.min().data(), box.max().data(),
        [&includes, this, myci](int i) {
        if (_segments[i][2] != myci)
          includes.push_back(i);
        return true;
      });
      if (includes.empty()) 
        continue;

      // group segments according to curve id.
      boost::unordered_map<int, DIntArray> includesmap;
      for (int ii : includes) {
        int curveid = _segments[ii][2];
        includesmap[curveid].push_back(ii);
      }

      for (auto& it : includesmap) {
        int curveid = it.first;
        for (int si : it.second) {
          // check parallelism, calculate intersection.
          Point3D pti0(_curvenodes[_segments[si][0]].data());
          Point3D pti1(_curvenodes[_segments[si][1]].data());
          Vector3D dir1(pti0, pti1);
          double len1 = dir1.length();
          if (len1 < Gloabals::ZERO) 
            continue;
          dir1 /= len1;
          if (std::abs(dir0.dot(dir1)) > angtol) 
            continue;
          Segment3D segi(pti0, pti1);
          double sc(-1.0), tc(-1.0);
          double dis2 = seg0.distance2(segi, sc, tc);
          if (dis2 > rad2) 
            continue;
          if (sc<0.001 || sc>(1.0 - 0.001)) 
            continue;
          // add to min dis.
          auto& mindisarray = mindismap[curveid];
          if (mindisarray.empty()) mindisarray.resize(nt - 1);
          mindisarray[ti] = std::min(mindisarray[ti], ItsctInfo(dis2, sc));
        }
      }
    }//for (int ti = 0; ti < nt - 1; ++ti)

  };

  auto lmdGetIntersections = [&]() {

    auto lmdAddIntersection = [&mytarray, &intersections, myci, &edge, this](int si, double sc, int destci) {
      double t = mytarray[si] + sc * (mytarray[si + 1] - mytarray[si]);
      Vector3D dir0, dir1;
      edge.tangent(t, dir0);
      Point3D pt0 = edge.eval(t);
      Point3D ptonedge;
      double tonedge;
      _curvedataarray[destci].edge.project(pt0.data(), ptonedge.data(), tonedge);
      _curvedataarray[destci].edge.tangent(tonedge, dir1.data());
      if (std::abs(dir0.dot(dir1)) < 0.985)
        intersections.push_back(t);
    };

    for (auto& mindisarrayit : mindismap) {
      double invalid0 = 0.01*Gloabals::MAX_REAL;
      double invalid1 = 0.1*Gloabals::MAX_REAL;
      double mindis = invalid1;
      int minsi = -1;
      int ns = mindisarrayit.second.size();
      for (int si = 0; si < ns; ++si) {
        double dis = mindisarrayit.second.at(si).dis;
        if (dis < mindis) {
          mindis = dis;
          minsi = si;
        }
        // intersection is found.
        if (dis > invalid0) {
          if (minsi > -1)
            lmdAddIntersection(minsi, mindisarrayit.second.at(minsi).sc, mindisarrayit.first);
          mindis = invalid1;
          minsi = -1;
        }
      }
    }// for (auto& mindisarrayit : mindismap)

  };

  lmdCollectMinDis();
  lmdGetIntersections();

}

bool AlgoPASM::is_valid_fixed_edge(const Int2& edge0, MyShellMesh& msh, double rCloseEdge)
{
  Segment2D seg0(msh.get_node_coord2D(edge0[0]), msh.get_node_coord2D(edge0[1]));
  Segment3D seg3d0(msh.get_node_coord3D(edge0[0]), msh.get_node_coord3D(edge0[1]));
  double len2 = msh.get_node_coord3D(edge0[0]).distance2(msh.get_node_coord3D(edge0[1]));
  IntSet checked;
  bool hasclose(false);
  std::function<bool(const Int2&)> IsIntersected = 
    [&msh, &seg0, &seg3d0, &edge0, &hasclose, &checked, rCloseEdge, len2, this](const Int2& lnk) {
    if (hasclose) 
      return false;// will stop searching.

    Int2 edge1;
    msh.get_elem_edge(lnk, edge1[0], edge1[1]);
    if (edge1[0] == edge0[1] || edge1[1] == edge0[1] ||
      edge1[0] == edge0[0] || edge1[1] == edge0[0])
      return false;

    int rt0 = msh.valid_tri(edge0[0], edge0[1], edge1[0]);
    int rt1 = msh.valid_tri(edge0[0], edge0[1], edge1[1]);
    if (rt0*rt1 > 0) 
      return false;

    rt0 = msh.valid_tri(edge1[0], edge1[1], edge0[0]);
    rt1 = msh.valid_tri(edge1[0], edge1[1], edge0[1]);

    if (rt0*rt1 > 0)
      return false;

    // intersected with existing fix edges.
    if (is_fixed_edge(msh, edge1[0], edge1[1])) {
      hasclose = true;
      return false;
    }

    // check node-edge proximity.
    Point3D pt3d;
    Point2D pt2d;
    for (int k = 0; k < 2; ++k) {
      if (msh.node_is(edge1[k], MyShellMesh::NODE_SUPER))
        continue;
      if (checked.find(edge1[k]) != checked.end())
        continue;
      checked.insert(edge1[k]);

      // check in 2D.
      msh.get_node_coord2D(edge1[k], pt2d);
      if (seg0.isOnSegment(pt2d, 0.001)) {
        hasclose = true;
        return false;
      }

      // check in 3D.
      msh.get_node_coord3D(edge1[k], pt3d);
      double dis2 = seg3d0.distance2(pt3d);
      double nodesize2 = msh.get_node_size(edge1[k])*msh.get_node_size(edge1[k]);
      nodesize2 = std::min(len2, nodesize2);
      if (dis2 < nodesize2*rCloseEdge*rCloseEdge) {
        hasclose = true;
        return false;
      }
    }
    return true;
  };

  Int2 fst(-1, -1);
  bool recoveryed(false);
  std::function<bool(const Int2&)> FindFirstTri = 
    [&msh, &edge0, &fst, &recoveryed, &IsIntersected](const Int2& lnk) {
    Int2 edge1;
    msh.get_elem_edge(lnk, edge1[0], edge1[1]);
    if (edge1[0] == edge0[1] || edge1[1] == edge0[1]) {
      recoveryed = true;
      return false;
    }
    if (IsIntersected(lnk)) {
      Int2 adj(-1, -1);
      bool rt = msh.get_elem_adj(lnk[0], lnk[1], adj);
      fst = adj;
      return false;
    }
    return true;
  };

  DInt2Array eq;
  fst.set(-1, -1);
  msh.visit_node_elem_ball(edge0[0], FindFirstTri);
  if (recoveryed)
    return true;
  if (hasclose)
    return false;
  if (fst[0] == -1)
    return false;
  eq.push_back(fst);
  msh.elem_broad_first_search(eq, IsIntersected);
  if (hasclose)
    return false;

  return true;
}

bool AlgoPASM::is_fixed_edge(const MyShellMesh&msh, int nd0, int nd1)
{
  if (!msh.node_is(nd0, MyShellMesh::NODE_EDGE) || !msh.node_is(nd1, MyShellMesh::NODE_EDGE))
    return false;
  return msh.is_fix_edge(nd0, nd1);
}

int AlgoPASM::insert_surface_node(const Point2D& pt2d,
  const Point3D& pt3d,
  double srcsize,
  MyShellMesh& msh,
  Mesher& mesher,
  std::function<void(int nd)> updateInfo,
  std::function<bool(int ele)> filter) 
{
  auto& locator = *(mesher.locator());
  Int2 whiTri;
  int loc = locator(pt2d.data(), whiTri, MyShellMesh::I0);
  if (Mesher::Locator::OUTSIDE == loc) {
    return -1;
  }
  if (Mesher::Locator::ONEDGE == loc) {
    int nd[2];
    msh.get_elem_edge(whiTri, nd[0], nd[1]);
    if (is_fixed_edge(msh, nd[0], nd[1])) {
      double tol2 = msh.distance2(nd[0], nd[1]);
      tol2 = tol2 * _alpha_nn*_alpha_nn;
      for (int k = 0; k < 2; ++k) {
        double dis2 = msh.get_node_coord3D(nd[k]).distance2(pt3d);
        if (dis2 < tol2)
          return nd[k];
      }
      return -1;
    }
  }
  else if (Mesher::Locator::ONVERTEX == loc) {
    return msh.get_elem_node(whiTri[0], whiTri[1]);
  }

  // check inside a elem.
  if (!filter(whiTri[0])) {
    return -1;
  }
  // proximity check for nodes.
  int closetnd = -1;
  double mindis2 = Gloabals::MAX_REAL;
  auto ele = msh.get_elem(whiTri[0]);
  for (int k = 0; k < 3; ++k) {
    int nd = ele[k];
    if (msh.node_is(nd, MyShellMesh::NODE_SUPER))
      continue;
    double ndsz = msh.get_node_size(nd);
    double tol2 = std::min(ndsz, srcsize)*_alpha_nn;
    tol2 = tol2 * tol2;
    double dis2 = msh.get_node_coord3D(nd).distance2(pt3d);
    if (dis2 < tol2) {
      if (dis2 < mindis2) {
        mindis2 = dis2;
        closetnd = nd;
      }
    }
  }
  if (closetnd > -1)
    return closetnd;
  if (msh.elem_is(whiTri[0], MyShellMesh::ELE_OUT))
    return -1;

  // proximity check for edges.
  int nd0 = -1, nd1 = -1;
  Point3D pt0, pt1;
  for (int i = 0; i < 3; ++i) {
    msh.get_elem_edge(whiTri[0], i, nd0, nd1);
    if (!is_fixed_edge(msh, nd0, nd1)) 
      continue;
    msh.get_node_coord3D(nd0, pt0);
    msh.get_node_coord3D(nd1, pt1);
    Segment3D seg(pt0, pt1);
    double dis2 = seg.distance2(pt3d);
    double tol2 = pt0.distance2(pt1);
    tol2 = tol2 * _alpha_ne*_alpha_ne;
    if (dis2 < tol2)
      return -1;
  }

  // do insertion.
  MyShellMesh::Point_Type apt;
  MyShellMesh::SetCoord2D(apt, pt2d);
  int nodeH = -1;
  msh.elem_checkout(whiTri[0], MyShellMesh::ELE_RIGID);// unfrozen and force insertion.
  mesher.insert(apt, nodeH, loc, whiTri, Mesher::WITH_CONSTRAIN, -1);
  if (nodeH >= MyShellMesh::I0) {
    updateInfo(nodeH);
  }
  return nodeH;
}

void AlgoPASM::setup_asm()
{
    _surfdataarray.reserve(_faceMesh.size());
    for (auto& it : _faceMesh) {
      _surfdataarray.push_back(std::make_shared<SurfData>(*it));
    }

    int nsurf = (int)(_surfdataarray.size());
#pragma omp parallel for
    for (int fi = 0; fi < nsurf; ++fi) {
      SurfData& surfdata = *_surfdataarray[fi];

      surfdata.meshdata.isSucceed = 0;
      double faceuvscal = surfdata.meshdata.uvScale;
      int fidx = _surfdataarray[fi]->face.getIndex();

      // generate background mesh.
      Int2IntMap tmp;
      if (MR_Enums::eSucceed != _PSLGFormFun(surfdata.eval, surfdata.bkgmesh, tmp))
        continue;
      if (MR_Enums::eSucceed != _bkgSurfMeshFun(surfdata.eval, surfdata.bkgmesh, tmp))
        continue;

      // prepare mesh.
      if (MR_Enums::eSucceed != _PSLGFormFun(surfdata.eval, surfdata.meshdata, surfdata.meshdata.nodeMap))
        continue;
      if (MR_Enums::eSucceed != surfdata.mesher.domesh(Mesher::DO_DELAUNAY | Mesher::DO_RECOVERY | Mesher::DO_MARK_OUT))
        continue;
      surfdata.sizeMeshPtr = std::make_shared<SizeTree>(surfdata.face, _dfMax, _dfMin);
      surfdata.sizeMeshPtr->updateSize(surfdata.meshdata);

      if (_d_m > Gloabals::ZERO) {
        for (const auto& it : surfdata.meshdata.nodeMap) {
          surfdata.mergedNodeMap[it.second] = std::abs(it.first);
        }
      }

      surfdata.meshdata.isSucceed = 1;
    }

    // collect bgm for hash.
    for (int fi = 0; fi < nsurf; ++fi) {
      if (_surfdataarray[fi]->meshdata.isSucceed == 0)
        continue;
      MyShellMesh& msh = _surfdataarray[fi]->bkgmesh;
      int nele = msh.elem_num() + MyShellMesh::I0;
      for (int ei = MyShellMesh::I0; ei < nele; ++ei) {
        if (msh.elem_is(ei, MyShellMesh::ELE_DEAD))
          continue;
        if (msh.elem_is(ei, MyShellMesh::ELE_OUT))
          continue;
        _triarray.push_back(Int2(fi, ei));
      }
    }

    // reset mark.
    for (int fi = 0; fi < nsurf; ++fi) {
      _surfdataarray[fi]->meshdata.isSucceed = 0;
    }

    if (_d_m > Gloabals::ZERO) {
      _mergedNodeId = _maxId;
    }

}

void AlgoPASM::get_node_projections(const Point3D& srcpt,
  double srcsize,
  bool checknormal,
  int mynd,
  int mysi,
  SurfNodeMap& surfInsertedNodeMap) 
{
  Box3Hash trihash(_triarray.size(), [this](size_t ti, Box3HashBox& box) {
    box.reset();
    Int2 atri = _triarray[ti];
    MyShellMesh& msh = _surfdataarray[atri[0]]->bkgmesh;
    Point3D pt;
    for (int i = 0; i < 3; ++i) {
      int nd = msh.get_elem_node(atri[1], i);
      msh.get_node_coord3D(nd, pt);
      box |= pt;
    }
    return true;
  });

  const double eps = 2.0*_d_s;
  const double eps2 = eps * eps;
  auto lmdGetInitialGuess = [&srcpt, eps2, checknormal, this](SurfData& surfdata,
    const Vector3D&nrml,
    const DIntArray&elems,
    Point3D& xyzio,
    Point2D& uvio) ->bool {

    double mindis2 = Gloabals::MAX_REAL;
    int minei = -1;
    Point3D pt3ds[3];
    Point3D closest, bary, minbary;
    for (int ei : elems) {
      auto ele = surfdata.bkgmesh.get_elem(ei);
      surfdata.bkgmesh.get_node_coord3D(ele[0], pt3ds[0]);
      surfdata.bkgmesh.get_node_coord3D(ele[1], pt3ds[1]);
      surfdata.bkgmesh.get_node_coord3D(ele[2], pt3ds[2]);
      Triangle3D tri(pt3ds[0], pt3ds[1], pt3ds[2]);
      double dis2 = tri.distance2(srcpt, closest, bary);
      if (dis2 < mindis2) {
        mindis2 = dis2;
        minei = ei;
        minbary = bary;
        xyzio = closest;
      }
    }

    if (mindis2 > eps2)
      return false;

    auto ele = surfdata.bkgmesh.get_elem(minei);
    Point2D pt2ds[3];
    surfdata.bkgmesh.get_node_coord2D(ele[0], pt2ds[0]);
    surfdata.bkgmesh.get_node_coord2D(ele[1], pt2ds[1]);
    surfdata.bkgmesh.get_node_coord2D(ele[2], pt2ds[2]);
    uvio = minbary[0] * pt2ds[0] + minbary[1] * pt2ds[1] + minbary[2] * pt2ds[2];

    return true;
  };

  auto lmdAddInsertedNodeMap = [this, &surfInsertedNodeMap, mynd](const Point3D& xyz, const Point2D& uv,
    double dis, int nd, int surfidx) {
    ProjInfo pdata;
    memcpy(pdata.data(), xyz.data(), 3 * sizeof(double));
    memcpy(pdata.data() + 3, uv.data(), 2 * sizeof(double));
    pdata[5] = dis;
    pdata[6] = 1.0*nd;
    {
      std::lock_guard<std::mutex> lock(_mutex);
      surfInsertedNodeMap[surfidx][mynd] = pdata;
    }
  };

  Vector3D nrml(0.0, 0.0, 0.0);
  auto lmdIsQualifiedProjection = [&srcpt, eps2, checknormal, this, &nrml](const Point3D& xyzio, double& dis2) {
    dis2 = srcpt.distance2(xyzio);
    if (dis2 > eps2)
      return false;
    if (checknormal && dis2 > _eps*_eps*100.0) {
      Vector3D vec(srcpt, xyzio);
      vec.normalize();
      if (std::abs(vec.dot(nrml)) < 0.866)// cos30.
        return false;
    }
    return true;
  };

  auto lmdMove2Inner = [&srcpt](SurfData& sdata, const Point2D& uvio0, const Point3D& xyzio0,
    Point2D& pt2d, Point3D& pt3d, int&nodeH) {
    auto& locator = *(sdata.mesher.locator());
    MyShellMesh& msh = sdata.meshdata;

    Int2 whiTri;
    int loc = locator(pt2d.data(), whiTri);
    if (Mesher::Locator::OUTSIDE == loc) {
      pt2d = uvio0;
      pt3d = xyzio0;
      return true;
    }
    else if (Mesher::Locator::ONVERTEX == loc) {
      nodeH = msh.get_elem_node(whiTri[0], whiTri[1]);
      return false;
    }

    if (!msh.elem_is(whiTri[0], MyShellMesh::ELE_OUT)) {
      return false;
    }

    // find the nearest node.
    double mindis2 = Gloabals::MAX_REAL;
    int mnd = -1;
    auto ele = msh.get_elem(whiTri[0]);
    for (int k = 0; k < 3; ++k) {
      int nd = ele[k];
      if (msh.node_is(nd, MyShellMesh::NODE_SUPER))
        continue;
      double dis2 = msh.get_node_coord3D(nd).distance2(pt3d);
      if (dis2 < mindis2) {
        mindis2 = dis2;
        mnd = nd;
      }
    }

    if (mnd < MyShellMesh::I0) {
      pt2d = uvio0;
      pt3d = xyzio0;
      return true;
    }
    else {
      if (mindis2 > xyzio0.distance2(srcpt)) {
        pt2d = uvio0;
        pt3d = xyzio0;
        return true;
      }
      else {
        msh.get_node_coord2D(mnd, pt2d);
        msh.get_node_coord3D(mnd, pt3d);
        nodeH = mnd;
        return true;
      }
    }

  };

  Box3HashBox box;
  box |= srcpt;
  box.thicknessIt(_d_s * 1.5);
  DIntArray includes;
  trihash.getBoxIncludes(box.min().data(), box.max().data(), [&includes, this, mysi](int i) {
    if (_triarray[i][0] == mysi) 
      return true;
    includes.push_back(i);
    return true;
  });

  if (includes.empty())
    return;

  boost::unordered_map<int, DIntArray> surfelems;
  for (int i : includes) {
    auto& elems = surfelems[_triarray[i][0]];
    elems.push_back(_triarray[i][1]);
  }

  SurfData& mysurfdata = *_surfdataarray[mysi];
  Point2D srcpt2d = mysurfdata.meshdata.get_node_coord2D(mynd);
  if (checknormal) {
    mysurfdata.eval.normal(srcpt2d[0], srcpt2d[1], nrml.data());
  }
  for (auto& it : surfelems) {
    int surfidx = it.first;
    SurfData& surfdata = *_surfdataarray[surfidx];
    TpFace& face = surfdata.face;

    Point3D xyzio0;
    Point2D uvio0(0.0, 0.0);
    if (!lmdGetInitialGuess(surfdata, nrml, it.second, xyzio0, uvio0))
      continue;

    Point3D xyzio(srcpt);
    Point2D uvio(uvio0);
    surfdata.eval.project(xyzio, uvio, TpFace::ePeojectByOptimize);

    double dis2 = -1.0;
    if (!lmdIsQualifiedProjection(xyzio, dis2))
      continue;

    int nodeH(0);
    if (lmdMove2Inner(surfdata, uvio0, xyzio0, uvio, xyzio, nodeH)) {
      if (!lmdIsQualifiedProjection(xyzio, dis2))
        continue;
      lmdAddInsertedNodeMap(xyzio, uvio, sqrt(dis2), nodeH, surfidx);
    }
    else {
      lmdAddInsertedNodeMap(xyzio, uvio, sqrt(dis2), 0, surfidx);
    }
  }

}

AlgoPASM::ImprintFailCode AlgoPASM::insert_edge(SurfData& srcdata,
  int mysi,
  const Int2& ei0,
  SurfData& dstdata,
  NodeProjectionMap& insertednodemap)
{
  MyShellMesh& srcmsh = srcdata.meshdata;

  auto lmdIsEdgeMatch = [this](SurfData& srcdata, const Int2& ei0, SurfData& dstdata, const Int2& ei1) {
    Point2D mid2d0 = 0.5*(srcdata.meshdata.get_node_coord2D(ei0[0]) + srcdata.meshdata.get_node_coord2D(ei0[1]));
    Point2D mid2d1 = 0.5*(dstdata.meshdata.get_node_coord2D(ei1[0]) + dstdata.meshdata.get_node_coord2D(ei1[1]));
    Point3D mid3d0, mid3d1;
    srcdata.eval(mid2d0, mid3d0);
    dstdata.eval(mid2d1, mid3d1);
    return mid3d0.distance2(mid3d1) < _d_s*_d_s*1.01;
  };

  auto lmdUpdataGId4Merge = [this](Int2IntMap& map0,
    Int2IntMap& map1,
    MyShellMesh& srcmsh,
    MyShellMesh& dstmsh,
    int ei0i, int ei1i) {
    std::lock_guard<std::mutex> lock(_mutex);

    auto gid0 = map0.find(ei0i);
    auto gid1 = map1.find(ei1i);
    if (gid0 == map0.end() && gid1 == map1.end()) {
      map0[ei0i] = ++_mergedNodeId;
      map1[ei1i] = _mergedNodeId;
    }
    else if (gid0 == map0.end() && gid1 != map1.end()) {
      map0[ei0i] = gid1->second;
    }
    else if (gid0 != map0.end() && gid1 == map1.end()) {
      map1[ei1i] = gid0->second;
    }
    else if (gid0 != map0.end() && gid1 != map1.end()) {
      if (gid0->second != gid1->second) {
        _gidMergePairs.push_back(Int2(gid0->second, gid1->second));
      }
    }
  };

  MRAssert(ei0[0] != ei0[1]);
  Int2 ei1(-1, -1);
  // insert two nodes even one of them failed.
  for (int k = 0; k < 2; ++k) {
    auto it = insertednodemap.find(ei0[k]);
    if (it == insertednodemap.end())
      continue;

    if (it->second[6] < -0.5) {// failed in insertion.
      continue;
    }
    else if (it->second[6] < 0.5) {// not inserted.
      double srcsize = srcmsh.get_node_size(ei0[k]);
      Point2D pt2d(it->second.data() + 3);
      int nodeH = insert_surface_node(pt2d, 
        Point3D(it->second.data()),
        srcsize, 
        dstdata.meshdata, 
        dstdata.mesher,
        [&dstdata, &pt2d, this, srcsize](int nodeH) {
        Point3D pt3d;
        dstdata.eval(pt2d, pt3d);
        dstdata.meshdata.set_node_coord3D(nodeH, pt3d);
        double nodesize = dstdata.sizeMeshPtr->sizeAt(pt2d[0], pt2d[1], true);
        nodesize = std::min(nodesize, srcsize);
        MRAssert(nodesize > 0.0);
        dstdata.meshdata.set_node_size(nodeH, nodesize);
        dstdata.meshdata.node_checkin(nodeH, MyShellMesh::NODE_FRONT);
      }, 
        [](int) {return true; 
      });
      it->second[6] = nodeH * 1.0;
      if (nodeH < 0)
        continue;
    }

    if (_d_m > Gloabals::ZERO) {
      if (it->second[5] < _d_m) {
        lmdUpdataGId4Merge(_surfdataarray[mysi]->mergedNodeMap,
          dstdata.mergedNodeMap,
          srcmsh, dstdata.meshdata,
          ei0[k], (int)it->second[6]);
      }
    }

    ei1[k] = (int)(it->second[6]);
  }

  if (ei1[0] < MyShellMesh::I0 || ei1[1] < MyShellMesh::I0)
    return ImprintFailCode::FailInsert;
  // add fix edges.
  if (ei1[0] == ei1[1])
    return ImprintFailCode::DegenerateEdge;
  MRAssert(ei1[0] >= MyShellMesh::I0 && ei1[1] >= MyShellMesh::I0);
  if (!lmdIsEdgeMatch(srcdata, ei0, dstdata, ei1))
    return ImprintFailCode::NonMatchEdge;
  if (!is_fixed_edge(dstdata.meshdata, ei1[0], ei1[1])) {
    if (!is_valid_fixed_edge(ei1, dstdata.meshdata, _alpha_ne))
      return ImprintFailCode::InvalidFixEdge;
    dstdata.meshdata.add_fix_edge(ei1[0], ei1[1], -1);
    if (MR_Enums::eSucceed != dstdata.mesher.recovery(ei1[0], ei1[1], [](int) {}))
      return ImprintFailCode::FailRecovery;
  }
  std::lock_guard<std::mutex> lock(_mutex);
  srcmsh.add_fix_edge(ei0[0], ei0[1], -1); // should also freeze src edges.

  return ImprintFailCode::Succeed;
}

void AlgoPASM::synchronize_surfaces_around_edges(SurfData& srcdata, const DInt2Array& edges, int mysi)
{
  MyShellMesh& srcmsh = srcdata.meshdata;
  DIntArray nodes;
  nodes.reserve(edges.size());
  for (auto& it : edges) {
    nodes.push_back(it[0]);
    nodes.push_back(it[1]);
  }
  Unique(nodes);

  // collect projections.
  SurfNodeMap surfInsertedNodeMap;
  int nnode = nodes.size();
#pragma omp parallel for
  for (int ni = 0; ni < nnode; ++ni) {
    int node = nodes[ni];
    Point3D coord;
    double size;
    srcmsh.get_node_coord3D(node, coord.data());
    size = srcmsh.get_node_size(node);
    bool checknormal = !srcmsh.node_is(node, MyShellMesh::NODE_RIDGE);
    get_node_projections(coord, size, checknormal, node, mysi, surfInsertedNodeMap);
  }
  if (surfInsertedNodeMap.size() < 1)
    return;

  std::vector<SurfNodeMap::iterator> its;
  int nsurf = surfInsertedNodeMap.size();
  its.reserve(nsurf);
  for (auto it = surfInsertedNodeMap.begin(); it != surfInsertedNodeMap.end(); ++it) {
    its.push_back(it);
  }

#pragma omp parallel for
  for (int iti = 0; iti < nsurf; ++iti) {
    const int dstFaceIdx = its[iti]->first;
    int nimprinted = 0;
    for (auto& ei : edges) {
      auto rt = insert_edge(srcdata, mysi, ei, *_surfdataarray[dstFaceIdx], its[iti]->second);
      if (ImprintFailCode::Succeed == rt) {
        ++nimprinted;
      }
    }

    if (nimprinted < 1)
      continue;

    MyShellMesh& dstmsh = _surfdataarray[dstFaceIdx]->meshdata;

    // frozen elems.
    int nele = dstmsh.elem_num() + MyShellMesh::I0;
    for (int ei = MyShellMesh::I0; ei < nele; ++ei) {
      if (dstmsh.elem_is(ei, MyShellMesh::ELE_DEAD))
        continue;
      if (dstmsh.elem_is(ei, MyShellMesh::ELE_OUT))
        continue;
      if (dstmsh.elem_is(ei, MyShellMesh::ELE_RIGID))
        continue;
      int i = 0;
      int nd0, nd1;
      for (; i < 3; ++i) {
        dstmsh.get_elem_edge(ei, i, nd0, nd1);
        if (!is_fixed_edge(dstmsh, nd0, nd1))
          break;
      }
      if (i == 3) {
        dstmsh.elem_checkin(ei, MyShellMesh::ELE_RIGID);
      }
    }

  }

}

void AlgoPASM::synchronize_surfaces_around_boundary_edges(const std::vector<std::pair<int, int> >& accessor)
{
  DInt2Array bedges;
  int nsurf = _surfdataarray.size();
  for (int ai = 0; ai < nsurf; ++ai) {
    const int si = accessor[ai].first;

    int surfid = _surfdataarray[si]->face.getIndex();
    if ((si % (int)ceil(1.0*nsurf / 10.0)) == 0)
      Messager::Msg(Messager::eMessage, "meshing surface patch %d [ %d / %d ]", surfid, ai, nsurf);

    MeshData& msh = _surfdataarray[si]->meshdata;

    bedges.clear();
    msh.getMeshEdge(bedges,
      [&msh](int ele) {
      return !msh.elem_is(ele, MyShellMesh::ELE_OUT);
    },
      MyShellMesh::eBoundary);

    bedges.erase(std::remove_if(bedges.begin(), bedges.end(),
      [&msh, this](const Int2& it) {
      if (msh.distance2(it[0], it[1]) < _eps*_eps*0.01)
        return true;
      return false;
    }), bedges.end());

    synchronize_surfaces_around_edges(*_surfdataarray[si], bedges, si);
  }

}

void AlgoPASM::synchronize_surfaces_around_inner_edges(const std::vector<std::pair<int, int> >& accessor)
{
  int nsurf = _surfdataarray.size();
  for (int ai = 0; ai < nsurf; ++ai) {
    const int si = accessor[ai].first;

    int surfid = _surfdataarray[si]->face.getIndex();
    if ((si % (int)ceil(1.0*nsurf / 10.0)) == 0)
      Messager::Msg(Messager::eMessage, "meshing surface patch %d [ %d / %d ]", surfid, ai, nsurf);

    SurfData& surfdata = *_surfdataarray[si];
    MeshData& emsh = surfdata.meshdata;
    MyShellMesh& fmsh = surfdata.fullfacet;

    Int2Set oldedgeset;
    int nedge = emsh.edge_num() + MyShellMesh::I0;
    for (int ei = MyShellMesh::I0; ei < nedge; ++ei) {
      Int2 edge(emsh.get_edge(ei));
      oldedgeset.insert(edge);
    }
    const int oldnodenum = emsh.node_num();

    _surfMeshFun(surfdata.eval, emsh);

    DInt2Array inneredges;
    emsh.getMeshEdge(inneredges,
      [&emsh](int ei) { 
      return !emsh.elem_is(ei, MyShellMesh::ELE_OUT); 
    }, MyShellMesh::eInternal);

    inneredges.erase(std::remove_if(inneredges.begin(), inneredges.end(), [&oldedgeset, &emsh, this](const Int2& it) {
      if (oldedgeset.find(it) != oldedgeset.end())
        return true;
      if (emsh.distance2(it[0], it[1]) < _eps*_eps*0.01)
        return true;
      return false;
    }), inneredges.end());

    synchronize_surfaces_around_edges(surfdata, inneredges, si);

    emsh.isSucceed = 1;
  }

}

bool AlgoPASM::contract_node_pair()
{
    auto lmdUnifyMergePair = [&]() {
      if (_gidMergePairs.empty()) 
        return;

      for (auto& it : _gidMergePairs) {
        if (it[0] > it[1])
          std::swap(it[0], it[1]);
      }

      DIntArray uniqueGidMap(_mergedNodeId + 1);
      std::iota(uniqueGidMap.begin(), uniqueGidMap.end(), 0);
      for (;;) {
        bool touched = false;
        for (auto& it : _gidMergePairs) {
          int& lval = uniqueGidMap[it[1]];
          int& sval = uniqueGidMap[it[0]];
          if (lval > sval) {
            touched = true;
            lval = sval;
          }
        }
        if (!touched) break;
      }

      const int nsurf = _surfdataarray.size();
#pragma omp parallel for
      for (int si = 0; si < nsurf; ++si) {
        for (auto& gid : _surfdataarray[si]->mergedNodeMap) {
          gid.second = uniqueGidMap[gid.second];
        }
      }

    };

    if (_mergedNodeId < 0)
      return true;

    for (auto& it : _surfdataarray) {
      if (it->meshdata.elem_num() < 1)
        it->mergedNodeMap.clear();
    }

    _mergedNodeId = -1;
    for (auto& it : _surfdataarray) {
      for (auto& nd : it->mergedNodeMap) {
        _mergedNodeId = std::max(_mergedNodeId, nd.second);
      }
    }

    lmdUnifyMergePair();

    std::vector<Point3D> ptsnew(_mergedNodeId + 1, Point3D(0.0, 0.0, 0.0));
    DIntArray ptcnt(_mergedNodeId + 1, 0);
    int nsurf = _surfdataarray.size();
    for (int si = 0; si < nsurf; ++si) {
      SurfData& surfdata = *_surfdataarray[si];
      if (surfdata.mergedNodeMap.empty()) 
        continue;
      for (auto& it : surfdata.mergedNodeMap) {
        Point3D pt3d = surfdata.meshdata.get_node_coord3D(it.first);
        ptcnt[it.second]++;
        ptsnew[it.second] = ptsnew[it.second] + pt3d;
      }
    }

#pragma omp parallel for
    for (int si = 0; si < nsurf; ++si) {
      SurfData& surfdata = *_surfdataarray[si];
      for (auto& it : surfdata.mergedNodeMap) {
        MRAssert(ptcnt[it.second] > 0);
        Point3D coordnew = ptsnew[it.second] / (double)ptcnt[it.second];
        surfdata.meshdata.set_node_coord3D(it.first, coordnew);
      }
    }

    return true;
}

void AlgoPASM::aligned_surface_meshing()
{
  setup_asm();

  // sort surfaces.
  std::vector<std::pair<int, int> > accessor;
  const int nsurf = _surfdataarray.size();
  accessor.reserve(nsurf);
  for (int si = 0; si < nsurf; ++si) {
    accessor.push_back(std::make_pair(si, _surfdataarray[si]->bkgmesh.elem_num()));
  }

  std::sort(accessor.begin(), accessor.end(),
    [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
    return a.second < b.second;
  });

  synchronize_surfaces_around_boundary_edges(accessor);
  synchronize_surfaces_around_inner_edges(accessor);

#pragma omp parallel for
  for (int si = 0; si < nsurf; ++si) {
    _surfdataarray[si]->mesher.domesh(Mesher::DELETE_OUT);
  }
}

