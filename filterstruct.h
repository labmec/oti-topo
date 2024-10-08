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

    REAL ComputeFiltereddcdxe(TPZVec<REAL> &dcdxf, TPZVec<REAL>& xf, TPZVec<REAL>& elvolvec, TPZCompMesh* cmesh) {

        // TPZGeoMesh* gmesh = cmesh->Reference();
        const REAL xe = xf[findex];
        const REAL myvol = elvolvec[findex];
        REAL sum = 0.;
        if (fneighIndexHf.size() == 0) {
            return dcdxf[findex];
        }
        for (const auto& neighIndexHf : fneighIndexHf) {
            int64_t neighIndex = neighIndexHf.first;
            REAL Hfneigh = neighIndexHf.second;
            if (Hfneigh < 0) DebugStop(); // what happened? It is defined as rmin - distcenters and should only be here if positeve
            TPZCompEl* cel = cmesh->Element(neighIndex);
            if(!cel) continue;
            const int64_t celindex = cel->Index(); // geoel and compel indexes are not the same
            REAL xneigh = xf[celindex];
            REAL dcdxneigh = dcdxf[celindex];
            REAL volneigh = elvolvec[celindex];
            sum += Hfneigh * xneigh * dcdxneigh * volneigh;
        }
        
        sum /= (fsumHf*xe*myvol);
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
        TPZCompMesh* cmesh = gmesh->Reference();
        TPZCompEl* cel = cmesh->Element(findex);
        TPZGeoEl* geo = cel->Reference();
        TPZGeoElSide geoside(geo);
        TPZManVector<REAL,3> rootcenter(3,0.);
        geoside.CenterX(rootcenter);
        TPZManVector<REAL,3> centerneigh(3,0.);
        int ncornernodes = geo->NCornerNodes();
        std::set<int64_t> tocheck, included, checked;
        tocheck.insert(findex);
        while (tocheck.size()) {
            TPZCompEl* celcheck = cmesh->Element(*tocheck.begin());
            TPZGeoEl* gel = celcheck->Reference();
            tocheck.erase(tocheck.begin());
            TPZGeoElSide gelside(gel);
            gelside.CenterX(centerneigh);
            REAL d = dist(centerneigh,rootcenter);
            REAL Hf = rmin - d;
            // std::cout << "Hf " << Hf << " d " << d << " rmin " << rmin << std::endl;
            if (Hf > 0) {
                included.insert(celcheck->Index());
                fneighIndexHf.push_back(std::make_pair(celcheck->Index(),Hf));
                fsumHf += Hf;
                for (int iside = 0; iside < gel->NCornerNodes(); iside++) {
                    TPZGeoElSide nodeside(gel,iside);
                    TPZGeoElSide neighbour = nodeside.Neighbour();
                    if (neighbour.Element()->HasSubElement()) {
                        TPZStack<TPZGeoEl*> subels;                        
                        // Get all the youngest children and insert them in the tocheck set
                        neighbour.Element()->YoungestChildren(subels);
                        for (auto subel : subels) {
                            if (subel->Dimension() == gmesh->Dimension()) {
                                int64_t neighindex = subel->Reference()->Index();
                                if (checked.find(neighindex) == checked.end() && included.find(neighindex) == included.end()) {
                                    tocheck.insert(neighindex);
                                }
                            }
                        }
                        continue;
                    }
                    while(neighbour != nodeside) {
                        if (neighbour.Element()->Dimension() == gmesh->Dimension() && neighbour.Element()->Reference()) {
                            int64_t neighindex = neighbour.Element()->Reference()->Index();
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
                checked.insert(gel->Reference()->Index());
            }
        }

    }
};