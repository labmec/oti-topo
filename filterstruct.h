#include <vector>
#include <utility>
#include "pzreal.h"
#include "pzvec_extras.h"

struct FilterStruct {    
    std::vector<std::pair<int, REAL>> fneighIndexHf;
    REAL fsumHf;
    int findex;

    FilterStruct(int index) : findex(index), fsumHf(0.) {
    }
    FilterStruct() : findex(-1), fsumHf(0.) {
    }

    REAL ComputeFiltereddcdxe(TPZVec<REAL> &dcdxf, TPZVec<REAL>& xf){

        REAL xe = xf[findex];
        REAL sum = 0.;
        for (const auto& neighIndexHf : fneighIndexHf) {
            int neighIndex = neighIndexHf.first;
            REAL Hfneigh = neighIndexHf.second;
            REAL xneigh = xf[neighIndex];
            REAL dcdxneigh = dcdxf[neighIndex];
            sum += Hfneigh * xneigh * dcdxneigh;
        }

        sum /= (fsumHf*xe);
        return sum;
    }

    void Print(std::ostream &out) {
        out << "Filter index " << findex << " sumHf " << fsumHf << " neighIndexHf " << " nneigh " << fneighIndexHf.size() << std::endl;
        for (const auto& neighIndexHf : fneighIndexHf) {
            out << neighIndexHf.first << " " << neighIndexHf.second << "\n";
        }
        out << std::endl;
    }

    void ComputeNeighIndexHf(TPZGeoMesh* gmesh, REAL rmin) {
        TPZGeoEl* geo = gmesh->ElementVec()[findex];
        TPZGeoElSide geoside(geo);
        TPZManVector<REAL,3> rootcenter(3,0.);
        geoside.CenterX(rootcenter);
        TPZManVector<REAL,3> centerneigh(3,0.);
        int ncornernodes = geo->NCornerNodes();
        std::set<int64_t> tocheck, included, checked;
        tocheck.insert(findex);
        while (tocheck.size()) {
            TPZGeoEl* gel = gmesh->ElementVec()[*tocheck.begin()];
            tocheck.erase(tocheck.begin());
            TPZGeoElSide gelside(gel);
            gelside.CenterX(centerneigh);
            REAL d = dist(centerneigh,rootcenter);
            REAL Hf = rmin - d;
            // std::cout << "Hf " << Hf << " d " << d << " rmin " << rmin << std::endl;
            if (Hf > 0) {
                included.insert(gel->Index());
                fneighIndexHf.push_back(std::make_pair(gel->Index(),Hf));
                fsumHf += Hf;
                for (int iside = 0; iside < gel->NCornerNodes(); iside++) {
                    TPZGeoElSide nodeside(gel,iside);
                    TPZGeoElSide neighbour = nodeside.Neighbour();
                    if (neighbour.Element()->HasSubElement()) {
                        continue;
                    }
                    while(neighbour != nodeside) {
                        if (neighbour.Element()->Dimension() == gmesh->Dimension()) {
                            int64_t neighindex = neighbour.Element()->Index();
                            // std::cout << "Verifying neighbour " << neighindex << std::endl;
                            if (checked.find(neighindex) == checked.end() && included.find(neighindex) == included.end()) {
                                // std::cout << "Inserting neighbour " << neighindex << std::endl;
                                tocheck.insert(neighindex);
                                
                            }
                        }
                        neighbour = neighbour.Neighbour();
                    }
                }                
            }
            else {
                checked.insert(gel->Index());
            }
        }

    }
};