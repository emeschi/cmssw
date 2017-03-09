#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalOneStepAlgo.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//
#include "DataFormats/CaloRecHit/interface/CaloID.h"


void HGCalOneStepAlgo::populate(const HGCRecHitCollection& hits){
  //loop over all hits and create the Hexel structure, skip energies below ecut
  for (unsigned int i=0;i<hits.size();++i) {
    const HGCRecHit& hgrh = hits[i];
    DetId detid = hgrh.detid();
    int layer = rhtools_.getLayerWithOffset(detid);
    float thickness = -9999.;
    unsigned thickIndex = -1;
    float sigmaNoise = -9999.;
    if(dependSensor){
      if( layer <= 40 ) {
	thickness = rhtools_.getSiThickness(detid);
	if( thickness>99. && thickness<101.) thickIndex=0;
	else if( thickness>199. && thickness<201. ) thickIndex=1;
	else if( thickness>299. && thickness<301. ) thickIndex=2;
	else assert( thickIndex>0 && "ERROR - silicon thickness has a nonsensical value" );
      }
      if( layer <= 40 ) sigmaNoise = 0.001 * fcPerEle * nonAgedNoises[thickIndex] * dEdXweights[layer] / (fcPerMip[thickIndex] * thicknessCorrection[thickIndex]);
      else if( layer <=52 ) sigmaNoise = 0.001 * noiseMip * dEdXweights[layer];
      if(hgrh.energy() < ecut*sigmaNoise) continue; //this sets the ZS threshold at ecut times the sigma noise for the sensor
    }
    if(!dependSensor && hgrh.energy() < ecut) continue; 


    layer += int(HGCalDetId(detid).zside()>0)*(maxlayer+1);
    
    // determine whether this is a half-hexagon
    bool isHalf = rhtools_.isHalfCell(detid);    
    const GlobalPoint position( std::move( rhtools_.getPosition( detid ) ) );
    //here's were the KDNode is passed its dims arguments - note that these are *copied* from the Hexel
    if(thickness<0.)thickness = 0.;
    unsigned int zside = unsigned(position.z()>0);
    float reducedz = float(rhtools_.getLayerWithOffset(detid))/2.;
    points[zside].pts.emplace_back(hgrh,detid,isHalf,sigmaNoise,thickness,&rhtools_);
    if(points[zside].pts.size()==0){
      minpos[zside][0] = position.x(); minpos[zside][1] = position.y(); minpos[zside][2] = reducedz;
      maxpos[zside][0] = position.x(); maxpos[zside][1] = position.y(); maxpos[zside][2] = reducedz; 
    }else{
      minpos[zside][0] = std::min((float)position.x(),minpos[zside][0]);
      minpos[zside][1] = std::min((float)position.y(),minpos[zside][1]);
      minpos[zside][2] = std::min((float)reducedz,minpos[zside][2]);
      maxpos[zside][0] = std::max((float)position.x(),maxpos[zside][0]);
      maxpos[zside][1] = std::max((float)position.y(),maxpos[zside][1]);
      maxpos[zside][1] = std::max((float)reducedz,maxpos[zside][2]);      
    }
  }

}
// Create a vector of Hexels associated to one cluster from a collection of HGCalRecHits - this can be used 
// directly to make the final cluster list - this method can be invoked multiple times for the same event 
// with different input (reset should be called between events)
void HGCalOneStepAlgo::makeClusters()
{
  //used for speedy search 

  //  std::vector<KDTree> hit_kdtree(2);

  //  std::vector<std::array<float,2> > minpos(2*(maxlayer+1),{ {0.0f,0.0f} }),maxpos(2*(maxlayer+1),{ {0.0f,0.0f} });

  //  std::vector<std::vector<Hexel> > points(2*(maxlayer+1)); //a vector of vectors of hexels, one for each layer
  //@@EM todo: the number of layers should be obtained programmatically - the range is 1-n instead of 0-n-1...


  if (verbosity < hgcal::pINFO)
    {
      std::cout << "-------------------------------------------------------------" << std::endl;
      std::cout << "HGC Imaging algorithm invoked for " << std::endl;
      std::cout << "delta_c " << delta_c << " kappa " << kappa;
      //      if( doSharing ) std::cout << " showerSigma " << std::sqrt(sigma2);
      std::cout << std::endl;
    }


  //assign all hits in each layer to a cluster core or halo
  for (unsigned int i = 0; i < 2; ++i) {

    
    my_kd_tree_t hit_kdtree(3,points[i],nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
    hit_kdtree.buildIndex();

    double maxdensity = calculateLocalDensity(points[i],hit_kdtree);
    // std::cout << "layer " << i << " max density " << maxdensity 
    // 	      << " total hits " << points[i].size() << std::endl;
    calculateDistanceToHigher(points[i]);
    findAndAssignClusters(points[i],hit_kdtree,maxdensity);
    //    std::cout << "found " << nclusters << " clusters" << std::endl;
  }
  //make the cluster vector
}

std::vector<reco::BasicCluster> HGCalOneStepAlgo::getClusters(bool doSharing){

  reco::CaloID caloID = reco::CaloID::DET_HGCAL_ENDCAP;
  std::vector< std::pair<DetId, float> > thisCluster;
  for (unsigned int i = 0; i < current_v.size(); ++i){
    double energy = 0;
    Point position = calculatePosition(current_v[i]);    
    std::vector< Hexel >::iterator it;
    for (it = current_v[i].begin(); it != current_v[i].end(); it++)
      {
	energy += (*it).isHalo ? 0. : (*it).weight;
	thisCluster.emplace_back(std::pair<DetId, float>((*it).detid,((*it).isHalo?0.:1.)));
      };
    if (verbosity < hgcal::pINFO)
      { 
	std::cout << "******** NEW CLUSTER (HGCIA) ********" << std::endl;
	std::cout << "Index          " << i                   << std::endl;
	std::cout << "No. of cells = " << current_v[i].size() << std::endl;
	std::cout << "     Energy     = " << energy << std::endl;
	std::cout << "     Phi        = " << position.phi() << std::endl;
	std::cout << "     Eta        = " << position.eta() << std::endl;
	std::cout << "*****************************" << std::endl;
      }
    clusters_v.push_back(reco::BasicCluster(energy, position, caloID, thisCluster, 
					    algoId));
    thisCluster.clear();
  }

  return clusters_v; 
}  

math::XYZPoint HGCalOneStepAlgo::calculatePosition(std::vector<Hexel> &v){
  float total_weight = 0.;
  float x = 0.;
  float y = 0.;
  float z = 0.;
  for (unsigned int i = 0; i < v.size(); i++){
    if(!v[i].isHalo){
      total_weight += v[i].weight;
      x += v[i].x*v[i].weight;
      y += v[i].y*v[i].weight;
      z += v[i].z*v[i].weight;
    }
  }
  
  return math::XYZPoint( x/total_weight, 
			 y/total_weight, 
			 z/total_weight );
} 

double HGCalOneStepAlgo::distance(const Hexel &pt1, const Hexel &pt2){
  return std::sqrt(distance2(pt1,pt2));
}

double HGCalOneStepAlgo::distance2(const Hexel &pt1, const Hexel &pt2){
  const double dx = pt1.x - pt2.x;
  const double dy = pt1.y - pt2.y;
  const double dz = pt1.lz - pt2.lz;
  return (dx*dx + dy*dy + dz*dz);
}


double HGCalOneStepAlgo::calculateLocalDensity(HexelCloud &nd, my_kd_tree_t &lp){
  double maxdensity = 0.;
  double search_radius = static_cast<double>(delta_c*delta_c);
  nanoflann::SearchParams params;
  for(unsigned int i = 0; i < nd.pts.size(); ++i){
    const double query_pt[3] = {nd.pts[i].x,nd.pts[i].y,nd.pts[i].lz};
    std::vector<std::pair<size_t,double> >  ret_matches; //index and distance of neighbors within radius
    const size_t nMatches = lp.radiusSearch(&query_pt[0],search_radius, ret_matches, params);
    if (verbosity < hgcal::pINFO)
      {
	std::cout << "hit " << i << " found " << nMatches << " nearest neighbors " << std::endl;
      }
    nd.pts[i].nNeighbors = nMatches;
    for(unsigned int k = 0; k < ret_matches.size(); k++){
      unsigned int j = ret_matches[k].first;
      if (verbosity < hgcal::pINFO)
	{
	  std::cout << "hit i,j " << i << "," << j 
		    << " distance by Hexel " 
		    << distance2(nd.pts[i],nd.pts[j]) 
		    << " distance by kdtree " 
		    <<  ret_matches[k].second
		    << std::endl;
	  std::cout << "delta x" << nd.pts[i].x-nd.pts[j].x << std::endl;
	  std::cout << "delta y" << nd.pts[i].y-nd.pts[j].y << std::endl;
	  std::cout << "delta lz" << nd.pts[i].lz-nd.pts[j].lz << std::endl;
	}
      
      nd.pts[i].rho += nd.pts[j].weight;
    }
    
    if(nd.pts[i].rho > maxdensity) maxdensity = nd.pts[i].rho;
    
    if (verbosity < hgcal::pINFO)
      {
	std::cout << "hit " << i << " local density is " << nd.pts[i].rho << std::endl;
      }
  }
  return maxdensity;
}

double HGCalOneStepAlgo::calculateDistanceToHigher(HexelCloud &nd){
  

  //sort vector of Hexels by decreasing local density
  std::vector<size_t> rs = sorted_indices(nd.pts);

  double maxdensity = 0.0;
  int nearestHigher = -1;


  if(rs.size()>0) 
    maxdensity = nd.pts[rs[0]].rho;
  else
    return maxdensity; // there are no hits
  double dist2 = 2500.0;
  //start by setting delta for the highest density hit to 
  //the most distant hit - this is a convention

  for(unsigned int j = 0; j < nd.pts.size(); j++){
    double tmp = distance2(nd.pts[rs[0]], nd.pts[j]);
    dist2 = tmp > dist2 ? tmp : dist2;
  }
  nd.pts[rs[0]].delta = std::sqrt(dist2);
  nd.pts[rs[0]].nearestHigher = nearestHigher;

  //now we save the largest distance as a starting point
  
  const double max_dist2 = dist2;
  
  for(unsigned int oi = 1; oi < rs.size(); ++oi){ // start from second-highest density
    dist2 = max_dist2;
    unsigned int i = rs[oi];
    // we only need to check up to oi since hits 
    // are ordered by decreasing density
    // and all points coming BEFORE oi are guaranteed to have higher rho 
    // and the ones AFTER to have lower rho
    for(unsigned int oj = 0; oj < oi; oj++){ 
      unsigned int j = rs[oj];
      double tmp = distance2(nd.pts[i], nd.pts[j]);
      if(tmp <= dist2){ //this "<=" instead of "<" addresses the (rare) case when there are only two hits
	dist2 = tmp;
	nearestHigher = j;
      }
    }
    nd.pts[i].delta = std::sqrt(dist2);
    nd.pts[i].nearestHigher = nearestHigher; //this uses the original unsorted hitlist 
    if (verbosity < hgcal::pINFO)
      {
	std::cout << "hit " << i << " distance to higher is " << nd.pts[i].delta << std::endl;
	std::cout << "hit " << i << " nearest higher is " << nd.pts[i].nearestHigher << std::endl;
      }

  }
  return maxdensity;
}

int HGCalOneStepAlgo::findAndAssignClusters(HexelCloud &nd,my_kd_tree_t &lp, double maxdensity){
  nanoflann::SearchParams params;
  params.sorted = true;
  double search_radius = static_cast<double>(delta_c*delta_c);

  //this is called once per side...
  //so when filling the cluster temporary vector of Hexels we resize each time by the number 
  //of clusters found. This is always equal to the number of cluster centers...

  unsigned int clusterIndex = 0;

  std::vector<size_t> rs = sorted_indices(nd.pts); // indices sorted by decreasing rho
  std::vector<size_t> ds = sort_by_delta(nd.pts); // sort in decreasing distance to higher


  for(unsigned int i =0; i < nd.pts.size(); ++i){

    //    std::cout << " delta " << lp[ds[i]].delta << " rho " << lp[ds[i]].rho << std::endl;
    if(nd.pts[ds[i]].delta < delta_c) break; // no more cluster centers to be looked at 
    if(dependSensor){

      //float rho_c = std::min(kappa*nd[ds[i]].data.sigmaNoise,maxdensity/10.);
      float rho_c = kappa*nd.pts[ds[i]].sigmaNoise;

      if(nd.pts[ds[i]].rho < rho_c ) continue; // set equal to kappa times noise threshold
    }
    else if(nd.pts[ds[i]].rho < maxdensity/kappa  /* || lp[ds[i]].rho<0.001*/)
      continue; 
    

    nd.pts[ds[i]].clusterIndex = clusterIndex;
    if (verbosity < hgcal::pINFO)
      {
	std::cout << "Adding new cluster with index " << clusterIndex+cluster_offset << std::endl;
	std::cout << "Cluster center is hit " << ds[i] << std::endl;
      }
    clusterIndex++;
  }

  //at this point clusterIndex is equal to the number of cluster centers - if it is zero we are 
  //done
  if(clusterIndex==0) return clusterIndex;

  //assign to clusters, using the nearestHigher set from previous step (always set except 
  // for top density hit that is skipped...by definition it is the ceneter of cluster 0
  for(unsigned int oi =1; oi < nd.pts.size(); ++oi){
    unsigned int i = rs[oi];
    int ci = nd.pts[i].clusterIndex;
    if(ci == -1){
      nd.pts[i].clusterIndex =  nd.pts[nd.pts[i].nearestHigher].clusterIndex;
    }
  }

  //make room in the temporary cluster vector for the additional clusterIndex clusters 
  // from this layer
  if (verbosity < hgcal::pINFO)
    { 
      std::cout << "resizing cluster vector by "<< clusterIndex << std::endl;
    }
  current_v.resize(cluster_offset+clusterIndex);

  //assign points closer than dc to other clusters to border region
  //and find critical border density
  std::vector<double> rho_b(clusterIndex,0.);

  //now loop on all hits again :( and check: if there are hits from another cluster within d_c -> flag as border hit

  for(unsigned int i = 0; i < nd.pts.size(); ++i){
    int ci = nd.pts[i].clusterIndex;
    bool flag_isolated = true;
    if(ci != -1){
      const double query_pt[3] = {nd.pts[i].x,nd.pts[i].y,nd.pts[i].lz};
      std::vector<std::pair<size_t,double> >  ret_matches; //index and distance of neighbors within radius
      const size_t nMatches = lp.radiusSearch(&query_pt[0],search_radius, ret_matches, params);
      if(nMatches!=nd.pts[i].nNeighbors){
	std::cout << "for hit " << i << " of cluster index "<< ci << " neighbors found differ " 
		  << nMatches << " - "
		  << nd.pts[i].nNeighbors
		  << std::endl;
      }
      // for(unsigned int k = 0; k < ret_matches.size(); k++){
      // 	unsigned int j = ret_matches[k].first;
      // 	//check if the hit is not within d_c of another cluster
      // 	if(nd.pts[j].clusterIndex!=-1){
      // 	  float dist = distance(nd.pts[j],nd.pts[i]);
      // 	  if(dist < delta_c && nd.pts[j].clusterIndex!=ci){
      // 	    if (verbosity < hgcal::pINFO)
      // 	      {
      // 		std::cout << "for hit " << i << " of cluster index "<< ci << " found cluster at " 
      // 			  << nd.pts[j].x << ","
      // 			  << nd.pts[j].y << ","
      // 			  << nd.pts[j].lz << " with cluster index " << nd.pts[j].clusterIndex << std::endl;
      // 	      }

      // 	     //in which case we assign it to the border
      // 	    nd.pts[i].isBorder = true;
      // 	    //check if this border hit has density larger than the current rho_b and update
      // 	    if(nd.pts[i].isBorder && rho_b[ci] < nd.pts[i].rho){
      // 	      if (verbosity < hgcal::pINFO)
      // 		{
      // 		  std::cout << "for hit " << i << " of cluster index "<< ci << " updating rho_b: before " <<  rho_b[ci] 
      // 			    << " after " << nd.pts[i].rho << std::endl;
      // 		}
      // 	      rho_b[ci] = nd.pts[i].rho;
      // 	    }
      // 	    break;
      // 	  }
      // 	  if(i!=j && nd.pts[j].clusterIndex==ci){
      // 	    //this is not an isolated hit
      // 	    if (verbosity < hgcal::pINFO)
      // 	      {
      // 		std::cout << "for hit " << i << " of cluster index "<< ci << " found at least one hit " 
      // 			  << nd.pts[j].x << ","
      // 			  << nd.pts[j].y << ","
      // 			  << nd.pts[j].lz << " with same cluster index " << nd.pts[j].clusterIndex 
      // 			  << std::endl;
      // 	      }

      // 	    flag_isolated = false;
      // 	  }
      // 	}
      // }

      // USE THE FOLLOWING CRITERION INSTEAD: the matches are sorted (closest point first)
      // if the closest NN is from another cluster, mark the hit as border and update the 
      // critical border density, otherwise keep hit as core (but still check if isolated)
      if(nMatches != 1){
	assert(i==ret_matches[0].first);
	for(unsigned int k = 1; k < ret_matches.size(); k++){
	  unsigned int j = ret_matches[k].first;
	  if(nd.pts[j].clusterIndex!=-1 && nd.pts[j].clusterIndex==ci){
	    flag_isolated = false;
	    break;
	  }
	  if(nd.pts[j].clusterIndex!=-1 && nd.pts[j].clusterIndex!=ci && nd.pts[j].rho>nd.pts[i].rho){
	    if (verbosity < hgcal::pINFO)
	      {
		std::cout << "for hit " << i << " of cluster index "<< ci << " found nearest hit at " 
			  << nd.pts[j].x << ","
			  << nd.pts[j].y << ","
			  << nd.pts[j].lz << " with cluster index " << nd.pts[j].clusterIndex 
			  << " rho[i] " << nd.pts[i].rho 
			  << " rho[j] " << nd.pts[j].rho << std::endl;
	      }
	    //in which case we assign it to the border
      	    nd.pts[i].isBorder = true;
	  }
	}
      }
      if(flag_isolated){
	nd.pts[i].isBorder = true; //the hit is more than delta_c from any of its brethren
      }
      //check if this border hit has density larger than the current rho_b and update
      if(nd.pts[i].isBorder && rho_b[ci] < nd.pts[i].rho){
	if (verbosity < hgcal::pINFO)
	  {
	    std::cout << "for hit " << i << " of cluster index "<< ci << " updating rho_b: before " <<  rho_b[ci] 
		      << " after " << nd.pts[i].rho << std::endl;
	  }
	rho_b[ci] = nd.pts[i].rho;
      }

    }	  

  }

  //flag points in cluster with density < rho_b as halo points, then fill the cluster vector 
  for(unsigned int i = 0; i < nd.pts.size(); ++i){
    int ci = nd.pts[i].clusterIndex;
    if(ci!=-1 && nd.pts[i].rho < rho_b[ci]){
      if (verbosity < hgcal::pINFO)
	{
	  std::cout << "for hit " << i << " of cluster index "<< ci << " rho_b is " <<  rho_b[ci]
		    << " hit density is " << nd.pts[i].rho << std::endl;
	}
      nd.pts[i].isHalo = true;
    }
    if(nd.pts[i].clusterIndex!=-1){ 
      current_v[ci+cluster_offset].push_back(nd.pts[i]);
      if (verbosity < hgcal::pINFO)
	{ 
	  std::cout << "Pushing hit " << i << " into cluster with index " << ci+cluster_offset 
		    << " is border " << nd.pts[i].isBorder 
		    << " is halo " << nd.pts[i].isHalo << std::endl;
	  std::cout << "Size now " << current_v[ci+cluster_offset].size() << std::endl;
	}
    }
  }

  //prepare the offset for the next layer if there is one
  if (verbosity < hgcal::pINFO)
    { 
      std::cout << "moving cluster offset by " << clusterIndex << std::endl;
    }
  cluster_offset += clusterIndex;
  return clusterIndex;
}

