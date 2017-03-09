#include <string>
#include <vector>
#include <set>
#include <iostream>

#include "RecoLocalCalo/HGCalRecAlgos/interface/KDTreeLinkerAlgoT.h"

class point{
public:
  point():x(0.),y(0.),z(0.){}
  point(float ix, float iy, float iz): x(ix),y(iy),z(iz){}
  float gx()const{return x;}
  float gy()const{return y;}
  float gz()const{return z;}

private:
  float x;
  float y;
  float z;
};
typedef KDTreeLinkerAlgo<point,3> KDTree;
typedef KDTreeNodeInfoT<point,3> KDNode;

int main(){
  
  std::vector<KDNode> points;
  KDTree hit_kdtree;
  KDTreeCube bounds(0.f,10.f,0.f,10.f,0.f,10.f);

  points.emplace_back(point(0.f,0.f,0.f),0.f,0.f,0.f);
  points.emplace_back(point(1.f,1.f,1.f),1.f,1.f,1.f);
  points.emplace_back(point(0.5f,1.f,0.5f),0.5f,1.f,0.5f);
  points.emplace_back(point(0.1f,0.f,0.f),0.1f,0.f,0.f);
  points.emplace_back(point(1.1f,1.f,1.f),1.1f,1.f,1.f);
  points.emplace_back(point(0.51f,1.f,0.5f),0.51f,1.f,0.5f);
  points.emplace_back(point(0.f,0.f,0.5f),0.f,0.f,0.5f);
  points.emplace_back(point(1.f,1.f,1.5f),1.f,1.f,1.5f);
  points.emplace_back(point(0.5f,1.f,0.55f),0.5f,1.f,0.55f);
  points.emplace_back(point(0.1f,0.f,0.5f),0.1f,0.f,0.5f);
  points.emplace_back(point(1.1f,1.f,1.5f),1.1f,1.f,1.5f);
  points.emplace_back(point(0.51f,1.f,0.55f),0.51f,1.f,0.55f);
  points.emplace_back(point(0.f,3.f,0.f),0.f,3.f,0.f);
  points.emplace_back(point(1.f,1.3f,1.f),1.f,1.3f,1.f);
  points.emplace_back(point(0.5f,1.3f,0.5f),0.5f,1.3f,0.5f);
  points.emplace_back(point(0.1f,0.3f,0.f),0.1f,0.3f,0.f);
  points.emplace_back(point(1.1f,1.3f,1.f),1.1f,1.3f,1.f);
  points.emplace_back(point(0.51f,1.3f,0.5f),0.51f,1.3f,0.5f);
  points.emplace_back(point(0.f,0.3f,0.5f),0.f,0.3f,0.5f);
  points.emplace_back(point(1.f,1.3f,1.5f),1.f,1.3f,1.5f);
  points.emplace_back(point(0.5f,1.3f,0.55f),0.5f,1.3f,0.55f);
  points.emplace_back(point(0.1f,0.3f,0.5f),0.1f,0.3f,0.5f);
  points.emplace_back(point(1.1f,1.3f,1.5f),1.1f,1.3f,1.5f);
  points.emplace_back(point(0.51f,1.3f,0.55f),0.51f,1.3f,0.55f);
  hit_kdtree.build(points,bounds);
  KDTreeCube search_box(-0.1f,.25f,-0.1f,.25f,-0.1f,.25f);
  
  std::vector<KDNode> found;
  hit_kdtree.search(search_box,found);
  std::cout << "here found " << found.size() << std::endl;
  for(unsigned int j = 0; j < found.size(); j++){
    std::cout << "found " << found[j].data.gx() << " " << found[j].data.gy() << " " << found[j].data.gz()
	      << std::endl;
  }
}
