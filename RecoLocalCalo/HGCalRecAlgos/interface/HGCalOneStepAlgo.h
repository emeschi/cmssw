#ifndef RecoLocalCalo_HGCalRecAlgos_HGCalOneStepAlgo_h
#define RecoLocalCalo_HGCalRecAlgos_HGCalOneStepAlgo_h

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"

// C/C++ headers
#include <string>
#include <vector>
#include <set>


//#include "KDTreeLinkerAlgoT.h"
#include "nanoflann.hpp"




class HGCalOneStepAlgo 
{

  
 public:
  
  
 HGCalOneStepAlgo() : delta_c(0.), kappa(1.), ecut(0.), cluster_offset(0),
		      sigma2(1.0),
		      algoId(reco::CaloCluster::undefined),
		      verbosity(hgcal::pERROR){
 }

  HGCalOneStepAlgo(float delta_c_in, double kappa_in, double ecut_in,
		   //		   const CaloSubdetectorTopology *thetopology_p,
		   reco::CaloCluster::AlgoId algoId_in,
                   bool dependSensor_in,
                   std::vector<double> dEdXweights_in,
                   std::vector<double> thicknessCorrection_in,
                   std::vector<double> fcPerMip_in,
                   double fcPerEle_in,
                   std::vector<double> nonAgedNoises_in,
                   double noiseMip_in,
		   hgcal::VerbosityLevel the_verbosity = hgcal::pERROR) : delta_c(delta_c_in), kappa(kappa_in), 
							    ecut(ecut_in),    
							    cluster_offset(0),
							    sigma2(1.0),
							    algoId(algoId_in),
                                                            dependSensor(dependSensor_in),
							    dEdXweights(dEdXweights_in),
                                                            thicknessCorrection(thicknessCorrection_in),
                                                            fcPerMip(fcPerMip_in),
                                                            fcPerEle(fcPerEle_in),
                                                            nonAgedNoises(nonAgedNoises_in),
                                                            noiseMip(noiseMip_in),
							    verbosity(the_verbosity),
							    points(2),
							    minpos(2,{ {0.0f,0.0f,0.0f} }),
							    maxpos(2,{ {0.0f,0.0f,0.0f} })
  {
  }

  HGCalOneStepAlgo(float delta_c_in, double kappa_in, double ecut_in,
		   double showerSigma, 
		   //		   const CaloSubdetectorTopology *thetopology_p,
		   reco::CaloCluster::AlgoId algoId_in,
                   bool dependSensor_in,
                   std::vector<double> dEdXweights_in,
                   std::vector<double> thicknessCorrection_in,
                   std::vector<double> fcPerMip_in,
                   double fcPerEle_in,
                   std::vector<double> nonAgedNoises_in,
                   double noiseMip_in,
		   hgcal::VerbosityLevel the_verbosity = hgcal::pERROR) : delta_c(delta_c_in), kappa(kappa_in), 
							    ecut(ecut_in),    
							    cluster_offset(0),
							    sigma2(std::pow(showerSigma,2.0)),
							    algoId(algoId_in),
                                                            dependSensor(dependSensor_in),
							    dEdXweights(dEdXweights_in),
                                                            thicknessCorrection(thicknessCorrection_in),
                                                            fcPerMip(fcPerMip_in),
                                                            fcPerEle(fcPerEle_in),
                                                            nonAgedNoises(nonAgedNoises_in),
                                                            noiseMip(noiseMip_in),
							    verbosity(the_verbosity),
							    points(2),
							    minpos(2,{ {0.0f,0.0f,0.0f} }),
							    maxpos(2,{ {0.0f,0.0f,0.0f} })
  {
  }

  virtual ~HGCalOneStepAlgo()
    {
    }

  void setVerbosity(hgcal::VerbosityLevel the_verbosity)
    {
      verbosity = the_verbosity;
    }

  void populate(const HGCRecHitCollection &hits);
  // this is the method that will start the clusterisation (it is possible to invoke this method more than once - but make sure it is with 
  // different hit collections (or else use reset)
  void makeClusters();
  // this is the method to get the cluster collection out 
  std::vector<reco::BasicCluster> getClusters(bool);
  // needed to switch between EE and HE with the same algorithm object (to get a single cluster collection)
  void getEventSetup(const edm::EventSetup& es){ rhtools_.getEventSetup(es); }
  // use this if you want to reuse the same cluster object but don't want to accumulate clusters (hardly useful?)
  void reset(){
    current_v.clear();
    clusters_v.clear();
    cluster_offset = 0;
    for( std::vector< HexelCloud >::iterator it = points.begin(); it != points.end(); it++)
      {
        // for( std::vector<KDNode>::iterator jt = it->begin(); jt != it->end(); jt++)
        //   delete jt->data;
        it->pts.clear();
      }
    for(unsigned int i = 0; i < minpos.size(); i++)
      {
	minpos[i][0]=0.;minpos[i][1]=0.;
	maxpos[i][0]=0.;maxpos[i][1]=0.;
      }
  }
  /// point in the space
  typedef math::XYZPoint Point;

 private: 
  
  //max number of layers
  static const unsigned int maxlayer = 52;

  // The two parameters used to identify clusters
  float delta_c;
  double kappa;

  // The hit energy cutoff
  double ecut;

  // The current offset into the temporary cluster structure
  unsigned int cluster_offset;

  // for energy sharing
  double sigma2; // transverse shower size

  // The vector of clusters
  std::vector<reco::BasicCluster> clusters_v;

  hgcal::RecHitTools rhtools_;

  // The algo id
  reco::CaloCluster::AlgoId algoId;

  // various parameters used for calculating the noise levels for a given sensor (and whether to use them)
  bool dependSensor;
  std::vector<double> dEdXweights;
  std::vector<double> thicknessCorrection;
  std::vector<double> fcPerMip;
  double fcPerEle;
  std::vector<double> nonAgedNoises;
  double noiseMip;

  // The verbosity level
  hgcal::VerbosityLevel verbosity;

  
  struct Hexel {

    double x;
    double y;
    double z;
    double lz;
    bool isHalfCell;
    double weight;
    double fraction;
    DetId detid;
    double rho;
    double delta;
    int nearestHigher;
    bool isBorder;
    bool isHalo;
    int clusterIndex;
    unsigned int nNeighbors;
    float sigmaNoise;
    float thickness;
    const hgcal::RecHitTools *tools;

    Hexel(const HGCRecHit &hit, DetId id_in, bool isHalf, float sigmaNoise_in, float thickness_in, const hgcal::RecHitTools *tools_in) : 
      x(0.),y(0.),z(0.),lz(0.),isHalfCell(isHalf),
      weight(0.), fraction(1.0), detid(id_in), rho(0.), delta(0.),
      nearestHigher(-1), isBorder(false), isHalo(false), 
      clusterIndex(-1), nNeighbors(0), sigmaNoise(sigmaNoise_in), thickness(thickness_in), 
      tools(tools_in)
    {
      const GlobalPoint position( std::move( tools->getPosition( detid ) ) );
      
      weight = hit.energy();
      x = position.x();
      y = position.y();
      z = position.z();
      lz = double(tools->getLayerWithOffset(detid))/2.;
    }
    Hexel() : 
      x(0.),y(0.),z(0.),lz(0.),isHalfCell(false),
      weight(0.), fraction(1.0), detid(), rho(0.), delta(0.),
      nearestHigher(-1), isBorder(false), isHalo(false), 
      clusterIndex(-1),
      nNeighbors(0),
      sigmaNoise(0.),
      thickness(0.),
      tools(0)
    {}
    bool operator > (const Hexel& rhs) const { 
      return (rho > rhs.rho); 
    }
    
  };

  struct HexelCloud
  {
    
    std::vector<Hexel>  pts;
    
    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }
    
    // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
    inline double kdtree_distance(const double *p1, const size_t idx_p2,size_t /*size*/) const
    {
      const double d0=p1[0]-pts[idx_p2].x;
      const double d1=p1[1]-pts[idx_p2].y;
      const double d2=p1[2]-pts[idx_p2].lz;
      return d0*d0+d1*d1+d2*d2;
    }
    
    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline double kdtree_get_pt(const size_t idx, int dim) const
    {
      if (dim==0) return pts[idx].x;
      else if (dim==1) return pts[idx].y;
      else return pts[idx].lz;
    }
    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
    
  };

  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, HexelCloud > ,
    HexelCloud,
    3 /* dim */
    > my_kd_tree_t;


  // A vector of vectors of KDNodes holding  Hexels in a cluster - to be used to build CaloClusters of DetIds
  std::vector< std::vector<Hexel> > current_v;

  std::vector<size_t> sort_by_delta(const std::vector<Hexel> &v){
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    sort(idx.begin(), idx.end(),
	 [&v](size_t i1, size_t i2) {return v[i1].delta > v[i2].delta;});
    return idx;
  }

  std::vector<HexelCloud> points; //a vector of vectors of hexels, one for each z side

  //bounding box
  std::vector<std::array<float,3> > minpos;
  std::vector<std::array<float,3> > maxpos;


  //these functions should be in a helper class.
  double distance2(const Hexel &pt1, const Hexel &pt2); //distance squared
  double distance(const Hexel &pt1, const Hexel &pt2); //2-d distance on the layer (x-y)
  double calculateLocalDensity(HexelCloud &, my_kd_tree_t &); //return max density
  double calculateDistanceToHigher(HexelCloud &);
  int findAndAssignClusters(HexelCloud &, my_kd_tree_t &, double);
  math::XYZPoint calculatePosition(std::vector<Hexel> &);


 };

#endif
