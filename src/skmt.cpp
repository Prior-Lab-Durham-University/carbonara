#include "ktlMoleculeRandom.h"
#include "skmt.h"

skmt::skmt(){};

template <typename T> int skmt::sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int skmt::signedTetraVolume(point&a,point&b,point&c,point&d){
    point ba = b.dif(a);
    point ca = c.dif(a);
    point da = d.dif(a);
    return sgn(ba.cross(ca).dotprod(da)/6.0);
}

bool skmt::intersectLineTriangle(point&q1, point&q2, point&p1, point&p2, point&p3){
    int s1 = signedTetraVolume(q1, p1, p2, p3);
    int s2 = signedTetraVolume(q2, p1, p2, p3);
    if(s1 != s2){
        int s3 = signedTetraVolume(q1, q2, p1, p2);
        int s4 = signedTetraVolume(q1, q2, p2, p3);
        int s5 = signedTetraVolume(q1, q2, p3, p1);
        if(s3 == s4 && s4 == s5){
            return true;
        }
    }
    return false;
}



std::vector<point> skmt::getSKMTCurve(ktlMolecule&molin){
    int i=0;
    std::vector<std::vector<point>> newcurve;
    std::vector<std::vector<point>> splitcurve;
    splitcurve = molin.getSubsecCoordinates(i);
    for(int subsec=0;subsec<splitcurve.size()-1;subsec++){
        if(splitcurve[subsec].size()>2){
            std::vector<bool> checks;
            for(int idx=1;idx<splitcurve[subsec].size();idx++){
                point p1 = splitcurve[subsec][0];
                point p2 = splitcurve[subsec][idx];
                point p3 = splitcurve[subsec+1][0];
                for(int subi=0;subi<splitcurve.size();subi++){
                    for(int subj=1;subj<splitcurve[subi].size();subj++){
                        point q0 = splitcurve[subi][subj-1];
                        point q1 = splitcurve[subi][subj];
                        checks.push_back(intersectLineTriangle(q0,q1,p1,p2,p3)); 
                    }
                }
            }
            bool noneAreTrue = std::none_of(checks.begin(), checks.end(), [](bool b) { return b; });
            if(noneAreTrue){
                newcurve.push_back({splitcurve[subsec][0]});
            }
            else{
                int idx=1;
                std::vector<bool> checks2;
                while(idx<splitcurve[subsec].size()){
                    point p1 = splitcurve[subsec][0];
                    point p2 = splitcurve[subsec][idx];
                    point p3 = splitcurve[subsec+1][0];
                    for(int subi=0;subi<splitcurve.size();subi++){
                        for(int subj=1;subj<splitcurve[subi].size();subj++){
                            point q0 = splitcurve[subi][subj-1];
                            point q1 = splitcurve[subi][subj];
                            checks2.push_back(intersectLineTriangle(q0,q1,p1,p2,p3));
                        }
                    }
                    bool noneAreTrue2 = std::none_of(checks2.begin(), checks2.end(), [](bool b) { return b; });
                    if(noneAreTrue2==true){
                        idx+=1;
                        checks2.clear();
                    }
                    else{
                        std::vector<point> vecPoints;
                        vecPoints.push_back(splitcurve[subsec][0]);
                        vecPoints.push_back(splitcurve[subsec][idx]);
                        newcurve.push_back(vecPoints);
                        idx = splitcurve[subsec].size();
                    }  
                }
            }
        }
        else{
            newcurve.push_back({splitcurve[subsec][0]});
        }
    }
    std::vector<point> lastPoints;
    lastPoints.push_back(splitcurve[splitcurve.size()-1][0]);
    lastPoints.push_back(splitcurve[splitcurve.size()-1][splitcurve[splitcurve.size()-1].size()-1]);
    newcurve.push_back(lastPoints);
    std::vector<point> skmtcurve;
    for(int i=0;i<newcurve.size();i++){
        for(int j=0;j<newcurve[i].size();j++){
            skmtcurve.push_back(newcurve[i][j]);
        }
    }
    return skmtcurve;   
}