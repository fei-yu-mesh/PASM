#pragma once


typedef std::function<void(TpEdgeWrapper&, const Real2&, DRealArray&)> CurveMeshFun;
typedef SurfMeshData0<6> MyShellMesh;
typedef std::function<int(SurfEvalFunctor& tpface, MyShellMesh& msh, Int2IntMap& node_map)> BackgroundSurfaceMeshFun;
typedef std::function<int(SurfEvalFunctor& tpface, MyShellMesh& msh, Int2IntMap& node_map)> PSLGFormFun;
typedef std::function<int(SurfEvalFunctor& tpface, MyShellMesh& msh)> SurfaceMeshFun;

class AlgoPASM
{
public:
  AlgoPASM(EntityId modelID, size_t, const int* faces);
  void set_CurveMeshFun(CurveMeshFun fun)
  {
    _curveMeshFun = fun;
  }
  void set_BackgroundSurfaceMeshFun(BackgroundSurfaceMeshFun fun)
  {
    _bkgSurfMeshFun = fun;
  }
  void set_PSLGFormFun(PSLGFormFun fun)
  {
    _PSLGFormFun = fun;
  }
  void set_SurfaceMeshFun(SurfaceMeshFun fun)
  {
    _surfMeshFun = fun;
  }
  int operator()(const char* ctrlStr);

private:
  typedef AlgoTriDelaunay2D<MyShellMesh> Mesher;
  typedef tBoxHash<double, 3> Box3Hash;
  typedef Box3Hash::Box Box3HashBox;

  struct EdgeData
  {
    TpEdge edge;
    std::vector<std::pair<int, double> > pts;
  };

  struct MeshData :public SurfMeshData0<6> {
    Int2IntMap nodeMap;// [global id, local id]
    TpFace     face;
    int        isSucceed;
    double     uvScale;
  };

  struct SurfData {
    SurfData(MeshData& md) :meshdata(md),
      face(md.face),
      eval(face),
      bkgmesh(),
      mesher(meshdata) {
      faceTp = face.type();
      meshdata.uvScale = face.get_uv_scale();
      meshdata.isSucceed = 0;
    }
    SurfData(const SurfData&) = delete;
    SurfData& operator = (const SurfData&) = delete;

  public:
    TpFace face;
    SurfEvalFunctor eval;
    int faceTp;
    MyShellMesh bkgmesh;
    std::shared_ptr<SizeTree> sizeMeshPtr;
    MeshData& meshdata;
    MyShellMesh fullfacet;
    Mesher mesher;
    Int2IntMap mergedNodeMap;
  };

  typedef boost::unordered_map<int, EdgeData> EdgeMeshMap;
  typedef std::shared_ptr<MeshData> MeshDataPtr;

  typedef Tuple0<double, 7> ProjInfo;
  typedef boost::unordered_map<int, ProjInfo> NodeProjectionMap;
  typedef boost::unordered_map<int, NodeProjectionMap> SurfNodeMap;

  struct CNode
  {
    CNode() :_id(0), _t(0.0), _df(-1.0) {}
    CNode(double t) :_id(0), _t(t), _df(-1.0) {}
    CNode(int id, double t) :_id(id), _t(t), _df(-1.0) {}

    bool operator < (const CNode& src) const
    {
      return this->_t < src._t;
    }

    int _id;
    double _t;
    double _df;
  };
  typedef std::vector<CNode> CNodeArray;

  struct CurveData
  {
    CurveData(TpEdge& e) :edge(e), isMeshed(false)
    {
      double intv[2];
      edge.interval(intv);
      TpVertex vx[2];
      edge.getVertex(vx);
      if (!edge.isForward())
        std::swap(vx[0], vx[1]);
      tarray1.push_back(CNode(vx[0].getIndex(), intv[0]));
      tarray1.push_back(CNode(vx[1].getIndex(), intv[1]));
    }

    TpEdge edge;
    DRealArray tarray4Itsct;
    DRealArray tarray0;
    CNodeArray tarray1;
    bool isMeshed;
  };

  struct NodeInfo
  {
    CNode cnd;
    Point3D pt3d;
    Vector3D tangent;
    std::function<bool(int)> filter;
    std::shared_ptr<IntSet> cset;
  };

  enum struct ImprintFailCode
  {
    Succeed = 0
    , FailProject
    , FailInsert
    , DegenerateEdge
    , InvalidFixEdge
    , NonMatchEdge
    , FailRecovery
  };

private:
  void setup(const char* ctrlStr);
  void setup_acm();
  void aligned_curve_meshing();
  void confirm_curve_mesh();
  void synchronize_curves_for_inner_nodes(const tBoxHash<double, 3>& segmenthash);
  void synchronize_curves_for_vertices(const tBoxHash<double, 3>& segmenthash);
  void synchronize_curves_around_nodes(const tBoxHash<double, 3>& segmenthash, std::vector<NodeInfo>& nodesinfo);
  void do_synchronize_curves_around_nodes(const tBoxHash<double, 3>& segmenthash, const CNodeArray& nodesnew, int ci);
  double get_curve_mesh_size(CurveData& curvedata, double t);
  bool insert_curve_node(CNode& tonedge, TpEdge& edge, CNodeArray& tarray, double rate, bool accurate);
  void synchronize_curves_around_node(const tBoxHash<double, 3>& segmenthash,
    const Point3D& pt,
    const Vector3D& dir,
    const CNode& cnd,
    std::shared_ptr<IntSet> ciset,
    std::function<bool(int)> filter,
    std::vector<NodeInfo>& nodesinfo2);
  void insert_touching_points(const tBoxHash<double, 3>& segmenthash);
  void get_touching_points(const tBoxHash<double, 3>& segmenthash,
    TpEdgeWrapper& edge,
    int myci,
    const DRealArray& mytarray,
    DRealArray& intersections);
  void aligned_surface_meshing();
  inline bool is_fixed_edge(const MyShellMesh&msh, int nd0, int nd1);
  bool is_valid_fixed_edge(const Int2& edge0, MyShellMesh& msh, double rCloseNode = 0.01);
  double get_curve_align_radius(double size)
  {
    double val = size * _d_c_relative + _d_c;
    double upbd = _dfMax * 0.2;
    if (val > upbd)
      val = upbd;
    if (val < _d_s)
      val = _d_s;
    return val;
  }
  int insert_surface_node(const Point2D&,
    const Point3D&,
    double srcsize,
    MyShellMesh& surfdata,
    Mesher& mesher,
    std::function<void(int nd)> updateInfo,
    std::function<bool(int ele)> filter = [](int) {return true; });
  void setup_asm();
  void get_node_projections(const Point3D& srcpt, double srcsize, bool checknormal, int mynd, int mysi, SurfNodeMap& surfInsertedNodeMap);
  AlgoPASM::ImprintFailCode insert_edge(SurfData& srcdata, int mysi, const Int2& ei0, SurfData& dstdata, NodeProjectionMap& insertednodemap);
  void synchronize_surfaces_around_edges(SurfData& srcdata, const DInt2Array& edges, int mysi);
  void synchronize_surfaces_around_boundary_edges(const std::vector<std::pair<int, int> >& accessor);
  void synchronize_surfaces_around_inner_edges(const std::vector<std::pair<int, int> >& accessor);
  bool contract_node_pair();

private:
  ICADModel*  _modePtr;
  double      _dfMax, _dfMin, _angle, _angleDet;
  int         _maxId, _lastError;
  mutable std::mutex  _mutex;

  double _gamma;
  double _eps;
  double _alpha_c, _alpha_nn, _alpha_ne;
  double _d_c_relative, _d_c, _d_s, _d_m;

  int _mergedNodeId;
  int _CNodeMaxId;
  static constexpr double _auxNodeDet = 0.1;
  static constexpr double _angle4background = 10.0;
  double _dfMaxBKG;
  double _dfMinBKG;

  std::vector<CurveData> _curvedataarray;
  std::vector<Real5> _curvenodes;// [x,y,z,size average len,size min len]
  std::vector<Int3> _segments;// [node0 id,node1 id,curve id]
  IntSet _validVertexNodeSet;

  std::vector<std::shared_ptr<SurfData> > _surfdataarray;
  DInt2Array _triarray;// [surf id,local elem id]
  DInt2Array _gidMergePairs;

  EdgeMeshMap _edgeMesh;// <edge topo handle, array of points <vertex id,t> >.
  EdgeMeshMap _edgeBKGMesh;// <edge topo handle, array of points <vertex id,t> >.
  std::vector<MeshDataPtr> _faceMesh;
  boost::unordered_map<int, MeshDataPtr> _faceMeshMap;

  BackgroundSurfaceMeshFun _bkgSurfMeshFun;
  SurfaceMeshFun _surfMeshFun;
  PSLGFormFun _PSLGFormFun;
  CurveMeshFun _curveMeshFun;

};

