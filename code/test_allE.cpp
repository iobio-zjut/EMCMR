#include "MC.h"
#include "GetVdwRadius.h"
#include <time.h>
#include "Rotbuilder.h"
#include "CaDensityx.h"
#include "randomx.h"
#include <unistd.h>

using namespace std;

/*
random_device rdTRNG;
mt19937 mtRNG;
std::uniform_real_distribution<double> fxdist;
std::uniform_int_distribution<int> ixdist; 

double randf0and1(){
    //The line below isn't necessary since the fdist distribution should always be between 0 and 1.
    //But it's there (for now) due to clarity.
    fxdist.param(std::uniform_real_distribution<double>(0.0, 1.0).param());
    return fxdist(mtRNG);
}
int randIntCustom(int a, int b){ 
    ixdist.param(std::uniform_int_distribution<int>(a, b).param());
    return ixdist(mtRNG);
}  */

std::string
name2eltxx( std::string line ) {
    std::string atmid = line.substr(12,4);
    while ( !atmid.empty() && atmid[0] == ' ' ) atmid = atmid.substr(1,atmid.size()-1);
    while ( !atmid.empty() && atmid[atmid.size()-1] == ' ' ) atmid = atmid.substr(0,atmid.size()-1);
    std::string resname = line.substr(17,3);
    while ( !resname.empty() && resname[0] == ' ' ) resname = resname.substr(1,resname.size()-1);
    while ( !resname.empty() && resname[resname.size()-1] == ' ' ) resname = resname.substr(0,resname.size()-1);

    std::string type;

    if ( line.substr(0,4) == "ATOM" ) {
        type = atmid.substr(0,2);
        if ( isdigit(type[0]) ) {
            // sometimes non-standard files have, e.g 11HH
            if ( !isdigit(type[1]) ) type = atmid.substr(1,1);
            else type = atmid.substr(2,1);
        } else if ( (line[12] == ' ' && type!="Zn" && type!="Fe" && type!="ZN" && type!="FE")
                || isdigit(type[1]) ) {
            type = atmid.substr(0,1);     // one-character element
        }

        if ( resname.substr(0,2) == "AS" || resname[0] == 'N' ) {
            if ( atmid == "AD1" ) type = "O";
            if ( atmid == "AD2" ) type = "N";
        }
        if ( resname.substr(0,3) == "HIS" || resname[0] == 'H' ) {
            if ( atmid == "AD1" || atmid == "AE2" ) type = "N";
            if ( atmid == "AE1" || atmid == "AD2" ) type = "C";
        }
        if ( resname.substr(0,2) == "GL" || resname[0] == 'Q' ) {
            if ( atmid == "AE1" ) type = "O";
            if ( atmid == "AE2" ) type = "N";
        }
        if ( atmid.substr(0,2) == "HH" ) { // ARG
            type = "H";
        }
        if ( atmid.substr(0,2) == "HD" || atmid.substr(0,2) == "HE" || atmid.substr(0,2) == "HG" ) {
            type = "H";
        }
    } else {
        if ( isalpha(atmid[0]) ) {
            if ( atmid.size() > 2 && (atmid[2] == '\0' || atmid[2] == ' ') ) {
                type = atmid.substr(0,2);
            } else if ( atmid[0] == 'A' ) { // alpha prefix
                type = atmid.substr(1, atmid.size() - 1);
            } else {
                type = atmid.substr(0,1);
            }
        } else if ( atmid[0] == ' ' ) {
            type = atmid.substr(1,1); // one char element
        } else {
            type = atmid.substr(1,2);
        }

        if ( atmid == resname ) {
            type = atmid;
            if ( type.size() == 2 ) type[1] = toupper(type[1]);
        } else if ( resname == "ADR" || resname == "COA" || resname == "FAD" ||
                resname == "GPG" || resname == "NAD" || resname == "NAL" ||
                resname == "NDP" || resname == "ABA" )  {
            if ( type.size() > 1 ) type = type.substr(0,1);
        } else if ( isdigit(type[0]) ) {
            type = type.substr(1,1);
        } else if ( type.size() > 1 && isdigit(type[1]) ) {
            type = type.substr(0,1);
        } else if ( type.size() > 1 && isalpha(type[1]) ) {
            if ( type[0] == 'O' && type[1] == 'H' ) {
                type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
            } else if ( islower(type[1]) ) {
                type[1] = toupper(type[1]);
            }
        }
    }
    return type;
}

vector<string> string_splitxx(string const & in, char splitchar)
{
    vector<string> parts;
    int i=0,j=0;
    while( j!= std::string::npos)
    {
        j=in.find(splitchar,i);
        std::string const part = in.substr(i,j-i);
        parts.push_back( part);
        i=j+1;
    }
    return parts;
}
string Int_to_String(int n)
{

    ostringstream stream;

    stream<<n; //n为int类型

    return stream.str();
}

void writePDBcoordss(std::string filename, poseCoords &atmlist)
{
    std::ifstream inpdb(filename.c_str());
//  cout<< filename.c_str()<<endl;
    std::string buf;
    vector<string> tmp_str;
    tmp_str = string_splitxx(filename,'.');
    string filenamex = tmp_str[0] + "_out.pdb";
    std::ofstream outpdb;
    outpdb.open(filenamex.c_str());
    int tmpx=0;

    while ( std::getline(inpdb, buf ) ) {
//        if ( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM" ) continue;
    	if(buf.length()<54) continue;
        if ( buf.substr(0,4) =="ATOM" && buf.substr(13,1)!="H") {

	        poseCoord atom_i;
	        string tmp_A,tmp_B;

	        tmp_A = buf.substr(0,30);
	        atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
	//        tmp_B = buf.substr(54,6);
	//        atom_i.B_ = atof(buf.substr(60,6).c_str());
	    //    cout<<"YYY"<<endl;


	//        atom_i.elt_ = name2elt( buf ); // horrible hacky logic mapping name->elt (could use PDB fields 76-77 if on by default)
	//        if ( atom_i.elt_ == "H" ) continue;

	        outpdb<<setw(30)<<tmp_A;
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[0];
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[1];
	        outpdb<<setw(8)<<setiosflags(ios::fixed)<<setprecision(3)<<atmlist[tmpx].x_[2];
	//        outpdb<<setw(6)<<tmp_B;
	//        outpdb<<setw(6)<<setiosflags(ios::fixed)<<setprecision(2)<<atmlist[tmpx].B_;
	        outpdb<<endl;
	        tmpx = tmpx + 1;        
	//        atmlist.push_back( atom_i );
    	}
    }
    outpdb.close();
}

// quick and dirty PDB read where we only care about heavyatom locations and atom ids
void
readPDBcoordss(std::string filename, poseCoords &atmlist) {
    std::ifstream inpdb(filename.c_str());
//  cout<< filename.c_str()<<endl;
    std::string buf;

    int nn=0;
    while ( !inpdb.eof() && std::getline(inpdb, buf) ) {
//    	cout<<nn << " ";
    	if(buf.length()<54) continue;
//    	if ( buf.substr(0,4) =="ATOM" && buf.substr(12,4) ==" CA " ) {
    	if ( buf.substr(0,4) =="ATOM" && buf.substr(13,1) !="H" && buf.substr(13,3) !="OXT") {
	        poseCoord atom_i;

	        atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
	        atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
	        atom_i.B_ = 0.0;// atof(buf.substr(60,6).c_str())
	        atom_i.elt_ = buf.substr(12,4);

//	        nn=nn+1;
//	        cout<< nn<<" ";
	        atmlist.push_back( atom_i );
    	}
    }
}

vector<float> Generate_random_point(float A,float B,float r)
{
	vector<float> coordx;
	float u=A+(B-A)*rand()/(RAND_MAX+1.0) ; // [a,b] random number,(rand()%(b-a+1))+a
	float v=A+(B-A)*rand()/(RAND_MAX+1.0);
	float theta = 2.0*MX_PI*u;
	float phi=acos(2*v-1.0);
	float x=r*sin(phi)*sin(theta);
	float y=r*cos(phi)*sin(theta);
	float z=r*cos(theta);
	coordx.push_back(x);
	coordx.push_back(y);
	coordx.push_back(z);

	return coordx;
}

float randomAB(float A,float B)
{
	return A+(B-A)*rand()/(RAND_MAX);
 //   return A+(B-A)*randf0and1();
}

float Distance_point(vector<float> A,vector<float> B)
{
	float dis=0.0;
	for(int i=0;i<A.size();i++)
	{
		dis = dis + (A[i]-B[i])*(A[i]-B[i]);
	}
	dis=sqrt(dis);
	return dis;
}
float Distance_pointy(vector<vector<float> > A,vector<vector<float> > B)
{

}
float Distance_pointx(vector<vector<float> > A,vector<vector<float> > B)
{
	float dis=0.0;
	float a1 = sqrt((A[0][0]-A[1][0])*(A[0][0]-A[1][0])+(A[0][1]-A[1][1])*(A[0][1]-A[1][1])+(A[0][2]-A[1][2])*(A[0][2]-A[1][2]));
	float a2 = sqrt((A[0][0]-A[2][0])*(A[0][0]-A[2][0])+(A[0][1]-A[2][1])*(A[0][1]-A[2][1])+(A[0][2]-A[2][2])*(A[0][2]-A[2][2]));
	float a3 = sqrt((A[0][0]-A[3][0])*(A[0][0]-A[3][0])+(A[0][1]-A[3][1])*(A[0][1]-A[3][1])+(A[0][2]-A[3][2])*(A[0][2]-A[3][2]));
	float b1 = sqrt((B[0][0]-B[1][0])*(B[0][0]-B[1][0])+(B[0][1]-B[1][1])*(B[0][1]-B[1][1])+(B[0][2]-B[1][2])*(B[0][2]-B[1][2]));
	float b2 = sqrt((B[0][0]-B[2][0])*(B[0][0]-B[2][0])+(B[0][1]-B[2][1])*(B[0][1]-B[2][1])+(B[0][2]-B[2][2])*(B[0][2]-B[2][2]));
	float b3 = sqrt((B[0][0]-B[3][0])*(B[0][0]-B[3][0])+(B[0][1]-B[3][1])*(B[0][1]-B[3][1])+(B[0][2]-B[3][2])*(B[0][2]-B[3][2]));

/*	for(int i=0;i<A.size();i++)
	{
		for(int j=0;j<A[i].size();j++)
		{
//			dis = dis + (A[i][j]-B[i][j])*(A[i][j]-B[i][j]);	
			dis = dis + (a1-b1)*(a1-b1)+		
		}
		
	} */
	dis = (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) + (a3-b3)*(a3-b3);
	dis=sqrt(dis);
	return dis;
} 


vector<float> calculateCoordinatesx(vector<float> refA,vector<float> refB,vector<float> refC,float L,float ang,float di)
{
	vector<float> AV(3,0.0),BV(3,0.0),CV(3,0.0);
	AV=refA;
	BV=refB;
	CV=refC;

	vector<float> CA(3,0.0),CB(3,0.0); 
	CA[0]=AV[0]-CV[0];
	CA[1]=AV[1]-CV[1];
	CA[2]=AV[2]-CV[2];
	CB[0]=BV[0]-CV[0];
	CB[1]=BV[1]-CV[1];
	CB[2]=BV[2]-CV[2];

	float AX = CA[0];
	float AY = CA[1];
	float AZ = CA[2];

	float BX = CB[0];
	float BY = CB[1];
	float BZ = CB[2];

	float A = (AY*BZ)-(AZ*BY);
	float B = (AZ*BX)-(AX*BZ);
	float G = (AX*BY)-(AY*BX);

	float F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang * (3.1415926/180.0));
//	float F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang);

	float constx = sqrt((pow(B*BZ-BY*G,2))*(-(F*F)*(A*A+B*B+G*G)+(B*B*(BX*BX+BZ*BZ) + A*A*(BY*BY+BZ*BZ)- (2*A*BX*BZ*G) + (BX*BX+ BY*BY)*G*G - (2*B*BY)*(A*BX+BZ*G))*L*L));
	float denom = (B*B)*(BX*BX+BZ*BZ)+ (A*A)*(BY*BY+BZ*BZ) - (2*A*BX*BZ*G) + (BX*BX+BY*BY)*(G*G) - (2*B*BY)*(A*BX+BZ*G);

	float X= ((B*B*BX*F)-(A*B*BY*F)+(F*G)*(-A*BZ+BX*G)+constx)/denom;

	float Y=0.0,Z=0.0;
    if((B==0 || BZ==0) && (BY==0 || G==0))
    {
        float const1=sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(L-X)*(L+X)));
        Y= ((-A*B*X)+const1)/(B*B+G*G);
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G));
    }
    else
    {
        Y= ((A*A*BY*F)*(B*BZ-BY*G)+ G*( -F*pow(B*BZ-BY*G,2) + BX*constx) - A*( B*B*BX*BZ*F- B*BX*BY*F*G + BZ*constx)) / ((B*BZ-BY*G)*denom);
        Z= ((A*A*BZ*F)*(B*BZ-BY*G) + (B*F)*pow(B*BZ-BY*G,2) + (A*BX*F*G)*(-B*BZ+BY*G) - B*BX*constx + A*BY*constx) / ((B*BZ-BY*G)*denom);
	}
	vector<float> D(3,0.0);
	D[0] = X + CV[0];
	D[1] = Y + CV[1];
	D[2] = Z + CV[2];

	float temp=0.0;

	temp=Points2Dihedral(AV,BV,CV,D) * (180.0/3.1415926);
//	cout<<"temp: "<<temp<<endl;

	di = di - temp;
	vector<float> Dx(3,0.0);
	Dx[0] = D[0] - BV[0];
	Dx[1] = D[1] - BV[1];
	Dx[2] = D[2] - BV[2];

	di = di * 3.1415926/180.0;
	GroupRotationt(CV,BV,di,Dx);
	D[0] = Dx[0] + BV[0];
	D[1] = Dx[1] + BV[1];
	D[2] = Dx[2] + BV[2];

	return D;
}

bool mcsuperpose(vector<poseCoord> &pointsBx,ElectronDensity theDensityMapx,int pnum)
{
    float KT=0.01;
    float new_E =0.0,old_E=0.0;
    vector<poseCoord> best_model;
    best_model = pointsBx;
    float best_E=10.0;
    ObjexxFCL::FArray3D< float > rhoC0x;
    ObjexxFCL::FArray3D< float > inv_rho_mask0x;
    rhoC0x.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
    inv_rho_mask0x.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
    // find initial position for A and B
    float init0max=-10000.0,init1max=-10000.0,init2max=-10000.0;
    float init0min=10000.0,init1min=10000.0,init2min=10000.0;       
    best_E=10.0;
    vector<float> coor_pdb;
    coor_pdb=vector<float>(3,0.0);   
    for(int j=0;j<pnum;j++)
    {        
        coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
        coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
        coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];       
    }
    coor_pdb[0] = coor_pdb[0]/(pnum);
    coor_pdb[1] = coor_pdb[1]/(pnum);
    coor_pdb[2] = coor_pdb[2]/(pnum);   
//    rhoC0x.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());
//    inv_rho_mask0x.dimension(theDensityMapx.density.u1() , theDensityMapx.density.u2() , theDensityMapx.density.u3());       
    for ( int t=0; t<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++t ) {
        rhoC0x[t]=0.0;
        inv_rho_mask0x[t]=1.0;
    }  
    vector<float> del_ijx,atm_jx;
    del_ijx=vector<float>(3,0.0);
    atm_jx=vector<float>(3,0.0);
    int tmp_i=0; 
    string elt_i; 
//    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
//    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;
    vector<vector<float>> atm_idx(pnum,vector<float>(3,0.0));   
    float bond0max=-10000.0,bond1max=-10000.0,bond2max=-10000.0;
    float bond0min=10000.0,bond1min=10000.0,bond2min=10000.0;       
    for(int i=0; i<pointsBx.size();i++)
    {
//        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
//        if(pointsBx[i].elt_ == "CA")
//        {
        if(pointsBx[i].elt_ == "CA")
        {           
//            cout<<"RRRR: "<<i;
            vector<float> cartX1;
            vector<float> fracX1;
            elt_i = pointsBx[i].elt_;
            elt_i = elt_i[0];
            OneGaussianScattering sig_j = get_A( elt_i );
            float k = sig_j.k( theDensityMapx.effectiveB );
            float C = sig_j.C( k );
            if ( C < 1e-6 ) continue;   

            cartX1 = pointsBx[i].x_;           
            MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

            // the location of atom in grid ?
            atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0] - theDensityMapx.origin[0] + 1) , (double)theDensityMapx.grid[0]);
            atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1] - theDensityMapx.origin[1] + 1) , (double)theDensityMapx.grid[1]);
            atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2] - theDensityMapx.origin[2] + 1) , (double)theDensityMapx.grid[2]);   
    //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
    //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
    //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
    //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
            for(int z=1;z<=theDensityMapx.density.u3();z++)
            {
                atm_jx[2] = z;
                del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMapx.grid[2];
                if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                for(int y=1;y<=theDensityMapx.density.u2();y++)
                {
                    atm_jx[1] = y;
                    del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMapx.grid[1];
                    // wrap-around??
                    if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                    if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                    del_ijx[0] = 0.0;
                    vector<float> frac_tmpy;
                    MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                    if(square_len(frac_tmpy)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                        
                    for(int x=1;x<=theDensityMapx.density.u1();x++)
                    {
                        atm_jx[0] = x;
                        del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMapx.grid[0];
                        // wrap-around??
                        if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                        if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                        vector<float> cart_del_ij2;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                        float d2 = square_len(cart_del_ij2);
                        if(d2 > (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;
                    
                        float atm = C*exp(-k*d2);
//                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                        float sigmoid_msk = exp( d2 - (theDensityMapx.ATOM_MASK)*(theDensityMapx.ATOM_MASK)  );
                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                    //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                        float inv_msk = 1/(1+sigmoid_msk);
                        rhoC0x(x,y,z) += atm;
                        inv_rho_mask0x(x,y,z) *= (1 - inv_msk);

                        if(x>bond0max) bond0max = x;
                        if(y>bond1max) bond1max = y;
                        if(z>bond2max) bond2max = z;
                        if(x<bond0min) bond0min = x;
                        if(y<bond1min) bond1min = y;
                        if(z<bond2min) bond2min = z;                                               

//                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
//                                mask2(x,y,z) = 1.0; // problem?
//                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
//                                    float atm = C*exp(-k*d2);
//                                    rhoC(x,y,z) += atm;
//                                }
//                            }
                    }
                }
            } 
//            cout<<" "<<tmp_i<<" "<<endl;           
            tmp_i = tmp_i + 1;
        }    
//        }
    }        
    float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
    float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
    float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
//    allpdbmass=0.0;
//    allmrcmass=0.0;
//    mass_pdb=vector<float>(3,0.0);
//    mass_mrc=vector<float>(3,0.0);
    for(int x=int(bond0min);x<=bond0max;x++)
    {
       for(int y=int(bond1min);y<=bond1max;y++)
        {
            for(int z=int(bond2min);z<=bond2max;z++)
            {
                clc_x2 = rhoC0x(x,y,z);
//                clc_x2 = tmp_rhc;
                obs_x2 = theDensityMapx.density(x,y,z);
//                obs_x2 = tmp_den;
                eps_x2 = 1-inv_rho_mask0x(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                sumCO_i2 += eps_x2*clc_x2*obs_x2;
                sumO_i2  += eps_x2*obs_x2;
                sumO2_i2 += eps_x2*obs_x2*obs_x2;
                sumC_i2  += eps_x2*clc_x2;
                sumC2_i2 += eps_x2*clc_x2*clc_x2;
                vol_i2   += eps_x2;
            }
        }
    } 
/*   for ( int x=0; x<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++x ) {
        clc_x2 = rhoC0x[x];
//            clc_x2 = tmp_rhc;
        obs_x2 = theDensityMapx.density[x];
//            obs_x2 = tmp_den;
        eps_x2 = 1-inv_rho_mask0x[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
//            eps_x2 = 1.0;

        // SMOOTHED
        sumCO_i2 += eps_x2*clc_x2*obs_x2;
        sumO_i2  += eps_x2*obs_x2;
        sumO2_i2 += eps_x2*obs_x2*obs_x2;
        sumC_i2  += eps_x2*clc_x2;
        sumC2_i2 += eps_x2*clc_x2*clc_x2;
        vol_i2   += eps_x2;
    }    */
    varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
    varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
//    cout<<"sumO2_i2: "<<sumO2_i2<<" "<<sumO_i2<<" "<<vol_i2<<endl;
//    cout<<"varC_i2: "<<varC_i2<<" "<<varO_i2<<endl;
    if ( varC_i2 == 0 || varO_i2 == 0 ) {
        CC_i2 = 0;
    } else {
        CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
    }   
    old_E=1.0-CC_i2;
    if(old_E<best_E) best_E = old_E;
//    cout<<"new_CC2: "<<new_CC2<<endl;
    for(int iii=0;iii<100;iii++) // 500
    {
//        cout<<"trdisx,anglex,KT: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;
//        float asin_theta=2.0*RandomDoubleX(0.0,1.0)-1.0;
        float asin_theta=2.0*randf0and1()-1.0;
        float acos_theta=sqrt(1.0-asin_theta*asin_theta);
        float apha=2.0*PI*randf0and1();
        float awx=acos_theta*cos(apha);
        float awy=acos_theta*sin(apha);
        float awz=asin_theta;
        // Translation Vector
        float t0=1.0;
    //    float t0=1.0;

    //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

        // Rotation matrix
        float anggg=90.0;
        float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
    //    if(new_CC2>0.95) t0=3.0;
    /*    if(acc_rate>40 && acc_rate<50)
        {
//            acc_rate=0;
        //    anggg=float(rand()%90)+120.0;
            anggg=float(randIntCustom(0,90))+90.0;
       //     anggg=float(randIntCustom(0,120));
            float tmp_rv=(2.0*randf0and1()-1.0);
            if(tmp_rv>=0.0)
            {
              angle_rotategg=anggg; // rotate angle  
            }
            else
            {
                angle_rotategg=(-1.0)*anggg; // rotate angle  
            }
        //    t0=float(rand()%5)+1.0;
        //    t0=float(randIntCustom(0,4))+1.0;
        //    t0=randf0and1();
            cout<<"rara: ";
        //   t0=abs(float(randf0and1())-0.5);
        //    angle_rotategg = (RandomDoubleX(0.0,1.0)*2.0-1.0)*anggg;
//                angle_rotategg = anggg;
//                t0=RandomDoubleX(0.0,3.0);
        }  */
/*        if(acc_rate>10 && acc_rate<15)
        {
            t0=randf0and1()+2.0;
        } */
    //    if(acc_rate>50) acc_rate=0; 
    /*    else if(acc_rate>50 && acc_rate<60)
        {
//            acc_rate=0;
        //    anggg=float(rand()%90)+120.0;
            anggg=float(randIntCustom(0,60))+120.0;
            float tmp_rv=(2.0*randf0and1()-1.0);
            if(tmp_rv>=0.0)
            {
              angle_rotategg=anggg; // rotate angle  
            }
            else
            {
                angle_rotategg=(-1.0)*anggg; // rotate angle  
            }                
        //    t0=float(rand()%5)+1.0;
            t0=abs(float(randf0and1())-0.5);  
        //    t0=float(randIntCustom(0,4))+1.0;
        //    t0=0.1;              
        }        */    
        float t1=(randf0and1()*2.0-1.0)*t0+0.0;
        float t2=(randf0and1()*2.0-1.0)*t0+0.0;
        float t3=(randf0and1()*2.0-1.0)*t0+0.0;  
    //    cout<<"ang,t0: "<<angle_rotategg<<" ("<<t1<<","<<t2<<","<<t3<<")"<<endl;          
        
        float asin=sin(angle_rotategg*(PI/180.0));
        float acos=cos(angle_rotategg*(PI/180.0));
        float u[3][3];
        u[0][0]=acos+awx*awx*(1.0-acos);
        u[0][1]=awx*awy*(1.0-acos)-awz*asin;
        u[0][2]=awx*awz*(1.0-acos)+awy*asin;
        u[1][0]=awx*awy*(1.0-acos)+awz*asin;
        u[1][1]=acos+awy*awy*(1.0-acos);
        u[1][2]=awy*awz*(1.0-acos)-awx*asin;
        u[2][0]=awx*awz*(1.0-acos)-awy*asin;
        u[2][1]=awy*awz*(1.0-acos)+awx*asin;
        u[2][2]=acos+awz*awz*(1.0-acos);            
        // Rotation points
//      int rand_point = rand()%(pnumx);
//      int rand_point = rand()%(pnum);
//      axyz[0]=pointsx[rand_point].ca[0];
//      axyz[1]=pointsx[rand_point].ca[1];
//      axyz[2]=pointsx[rand_point].ca[2];
    //    int rand_point = rand()%6+83 ;
    //    int rand_point = randop ;
    //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
    //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
    //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
        vector<float> axyz;
        axyz=vector<float>(3,0.0);        
        axyz[0]= coor_pdb[0];  
        axyz[1]= coor_pdb[1];
        axyz[2]= coor_pdb[2];

    //    beg_p=rand_point; 
//    cout<<"XXX2"<<endl;  

        int tmp3=0;
        vector<poseCoord > fin_mat;
    //    fin_mat = tmp_mat;
    //    if(beg_p==0)
    //    {
        tmp3=0;
        vector<poseCoord>().swap(fin_mat);
        for(int j=0;j<pnum;j++)
        {
            fin_mat.push_back(pointsBx[3*j+0]);
            fin_mat.push_back(pointsBx[3*j+1]);
            fin_mat.push_back(pointsBx[3*j+2]);
            tmp3=tmp3+1;
        } 
    //    old_E = new_E ;// + old_clash + 0.1*d1/d2;// + olddis_p_m;
    //    old_dE = old_clash + olddis_p_m;                      
        vector<poseCoord > fin_matx;
        fin_matx=fin_mat;
        tmp3=0;
        for(int j=0;j<fin_mat.size();j++) // not last N
        {
            fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
            fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
            fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
            tmp3=tmp3+1;
        }
        vector<poseCoord>().swap(fin_matx);

        for ( int t=0; t<theDensityMapx.density.u1()*theDensityMapx.density.u2()*theDensityMapx.density.u3(); ++t ) {
            rhoC0x[t]=0.0;
            inv_rho_mask0x[t]=1.0;
        }  
        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0);
        tmp_i=0; 
    //    string elt_i; 
    //    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
    //    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0; 
        atm_idx=vector<vector<float>>(pnum,vector<float>(3,0.0));  
        bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
        bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;       
        for(int i=0; i<fin_mat.size();i++)
        {
    //        cout<<"YYYY: "<<pointsBx[i].elt_<<endl;;
    //        if(pointsBx[i].elt_ == "CA")
    //        {
            if(fin_mat[i].elt_ == "CA")
            {           
    //            cout<<"RRRR: "<<i;
                vector<float> cartX1;
                vector<float> fracX1;
                elt_i = fin_mat[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                float k = sig_j.k( theDensityMapx.effectiveB );
                float C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_mat[i].x_;           
                MatrixTimesTransVector(theDensityMapx.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMapx.grid[0] - theDensityMapx.origin[0] + 1) , (double)theDensityMapx.grid[0]);
                atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMapx.grid[1] - theDensityMapx.origin[1] + 1) , (double)theDensityMapx.grid[1]);
                atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMapx.grid[2] - theDensityMapx.origin[2] + 1) , (double)theDensityMapx.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMapx.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMapx.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMapx.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMapx.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMapx.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMapx.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMapx.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK)*(theDensityMapx.ATOM_MASK_PADDING+theDensityMapx.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMapx.ATOM_MASK)*(theDensityMapx.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC0x(x,y,z) += atm;
                            inv_rho_mask0x(x,y,z) *= (1 - inv_msk);

                            if(x>bond0max) bond0max = x;
                            if(y>bond1max) bond1max = y;
                            if(z>bond2max) bond2max = z;
                            if(x<bond0min) bond0min = x;
                            if(y<bond1min) bond1min = y;
                            if(z<bond2min) bond2min = z;                                               

    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                } 
    //            cout<<" "<<tmp_i<<" "<<endl;           
                tmp_i = tmp_i + 1;
            }    
    //        }
        } 
        sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
        sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
        clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
        float allpdbmass=0.0;
        float allmrcmass=0.0;
        vector<float> mass_pdb;
        mass_pdb=vector<float>(3,0.0);
        vector<float> mass_mrc;
        mass_mrc=vector<float>(3,0.0);
        for(int x=int(bond0min);x<=bond0max;x++)
        {
           for(int y=int(bond1min);y<=bond1max;y++)
            {
                for(int z=int(bond2min);z<=bond2max;z++)
                {
                    clc_x2 = rhoC0x(x,y,z);
    //                clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMapx.density(x,y,z);
    //                obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask0x(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                }
            }
        }

        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
        if ( varC_i2 == 0 || varO_i2 == 0 ) {
            CC_i2 = 0;
        } else {
            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
        }
        new_E=1.0-CC_i2; 
//        new_E = new_CC1 ;                                   
//        cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<<old_E<<" "<<new_E<<endl;
        float dE = new_E-old_E;
        if(new_E<old_E)
        {
            tmp3=0;
            coor_pdb=vector<float>(3,0.0);
            for(int j=0;j<pnum;j++)
            {
                pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                    
                tmp3=tmp3+1;
            }
            coor_pdb[0] = coor_pdb[0]/(pnum);
            coor_pdb[1] = coor_pdb[1]/(pnum);
            coor_pdb[2] = coor_pdb[2]/(pnum); 
            if(new_E<best_E)
            {
                best_E = new_E;
                vector<poseCoord>().swap(best_model);
                best_model = fin_mat;
            }
        //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;      
            old_E = new_E;            
        //    new_CC2 = new_CC1;
             
//            cout<<"old_CC1,new_CC1: "<<old_E<<" "<<new_E<<endl;         
        }else
        {
        //    float tmpx=rand()/double(RAND_MAX);
            float tmpx=randf0and1();
            float mc_v = exp(-dE/(KT));
            if(tmpx< mc_v) // CC must be >0
            {
                tmp3=0;
                coor_pdb=vector<float>(3,0.0);
                for(int j=0;j<pnum;j++)
                {
                    pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                    pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                    pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                    coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                    coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                    coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                         
                    tmp3=tmp3+1;
                } 
                coor_pdb[0] = coor_pdb[0]/(pnum);
                coor_pdb[1] = coor_pdb[1]/(pnum);
                coor_pdb[2] = coor_pdb[2]/(pnum); 
            //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                     
            //    new_CC2 = new_CC1;
                old_E = new_E;
//                cout<<"old_CC1,new_CC1:XX "<<old_E<<" "<<new_E<<endl;                    
            }                
        }           
    } 
    for(int j=0;j<pnum;j++)
    {
        pointsBx[3*j+0] = best_model[3*j+0];
        pointsBx[3*j+1] = best_model[3*j+1];
        pointsBx[3*j+2] = best_model[3*j+2];
    //    bb[j]=bby[tmp_tmx];
    }     
    bool mcsp=true;
    return mcsp;   
}

Model REMC_sampley(Model posex, string inputM,float MRC_R,float mapsp,char bindir[])
{
    // read Density map
    ElectronDensity theDensityMap;
    theDensityMap.preso= MRC_R;
    theDensityMap.pATOM_MASK =7.0;
    theDensityMap.pCA_MASK =7.0;
    theDensityMap.pforce_apix_on_map_load_= 0.0;
    theDensityMap.pWINDOW_ =1;
    theDensityMap.pscore_window_context_ = false ;
    theDensityMap.premap_symm_ = false ;
    theDensityMap.pde_edensity = false;
    theDensityMap.pnkbins_ = 0 ;  
    theDensityMap.readMRCandResize(inputM,MRC_R,mapsp);
    float mapspx=3.0;
    ElectronDensity theDensityMapx;
    theDensityMapx.preso= MRC_R;
    theDensityMapx.pATOM_MASK =7.0;
    theDensityMapx.pCA_MASK =7.0;
    theDensityMapx.pforce_apix_on_map_load_= 0.0;
    theDensityMapx.pWINDOW_ =1;
    theDensityMapx.pscore_window_context_ = false ;
    theDensityMapx.premap_symm_ = false ;
    theDensityMapx.pde_edensity = false;
    theDensityMapx.pnkbins_ = 0 ;      
    theDensityMapx.readMRCandResize(inputM,MRC_R,mapspx);

    // model parameter
    Model fin_pose;
    vector<Rot> pointsx;
    vector<poseCoord> pointsBx;
    vector<poseCoord> pointsBx0;
    vector<string> pointsrenm;

    // residue number
    int pnum =0,pnumx=0;   
    for(int i=0;i<posex.chains.size();i++)
    {
        Chain Chanx = posex.chains[i];
        for(int j=0;j<Chanx.residues.size();j++)
        {
            Residue Resdx = Chanx.residues[j];
            Rot rotx;
            pointsrenm.push_back(Resdx.resname);
//            cout<<"atom: "<<Resdx.resname<<endl;
            vector<float> sc_coord(3,0.0);
            vector<float> co_coord(3,0.0);
            int sc_num=0;
            for(int t=0;t<Resdx.atoms.size();t++)
            {
                poseCoord pcdx;
                Atom Atmx=Resdx.atoms[t];
                if(Atmx.atona != " CA " && Atmx.atona != " C  " && Atmx.atona != " N  " && Atmx.atona != " O  ")
                {
                    rotx.res_atm.push_back(Atmx.xyzVect);
                    sc_coord[0] = sc_coord[0] + Atmx.xyzVect[0];
                    sc_coord[1] = sc_coord[1] + Atmx.xyzVect[1];
                    sc_coord[2] = sc_coord[2] + Atmx.xyzVect[2];
                    vector<string> atmnmx = string_splitxx(Atmx.atona,' ');
                    string atmnmy = atmnmx[0];  // error
                    rotx.resatmnm.push_back(atmnmy);
                    sc_num = sc_num + 1;
                }               
                if(Atmx.atona == " CA ") {rotx.ca = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);pnumx = pnumx +1;}
                if(Atmx.atona == " C  ") {rotx.cc = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);}//pnumx = pnumx +1;
                if(Atmx.atona == " N  ") {rotx.cn = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);}//pnumx = pnumx +1;
                if(Atmx.atona == " O  ") {co_coord[0] = Atmx.xyzVect[0];co_coord[1] = Atmx.xyzVect[1];co_coord[2] = Atmx.xyzVect[2];}  

            }
            sc_coord[0] = sc_coord[0]/float(sc_num);
            sc_coord[1] = sc_coord[1]/float(sc_num);
            sc_coord[2] = sc_coord[2]/float(sc_num);
            pointsx.push_back(rotx);
            pointsBx[3*j+2].x_o = co_coord;
            pointsBx[3*j+1].x_sc = sc_coord;            
//          pointsy.push_back(roty);
            pnum = pnum +1;                                     
        }
    }
    cout<<"pnum,pnumx: "<<pnum<<" "<<pnumx<<endl;
    pointsBx0 = pointsBx;

    cout<<"binddie: "<<bindir<<endl;
    char oneline[600];
    sprintf(oneline,"%s/bnramaprob",bindir);
    vector<vector<vector<float>>> rlogduke;
    vector<vector<vector<int>>> ramaduke;    
    loadramadukebn(oneline,rlogduke,ramaduke); // torsion angle distribution

    char onelinex[600];
    sprintf(onelinex,"%s/newsgdistriaanew72.txt",bindir);
    vector<vector<vector<vector<float>>>> sgpos2;    
    loadsgpos2(onelinex,72,sgpos2);

    char oneliney[600];
    sprintf(oneliney,"%s/phipsiaass.txt",bindir);
    float *phipsiprob[8][20];  
    loadphipsiprob(oneliney,phipsiprob);  

    vector<boneinfo> bb(pnum);
    vector<poseCoord> proseq;
    proseq=pointsBx;
    int numbb=pnum;
    int tmp_tt=0;
    for(int j=0;j<pnum;j++)
    {
        bb[tmp_tt].indn = 3*tmp_tt+0;
        bb[tmp_tt].indca = 3*tmp_tt+1;
        bb[tmp_tt].indc = 3*tmp_tt+2; 
        string resng = pointsrenm[tmp_tt];
        char resny = resng[0];
        bb[tmp_tt].resid = aminoid(resny);        
        tmp_tt=tmp_tt+1;
    }     
    cout<<"xxExxx"<<endl;
    calcssennhoc(bb,pnum,pointsBx);    

    vector<float> axyz;
    axyz=vector<float>(3,0.0);
    float ang; 
    float angle_rotate; 
    float old_dE1=0.0;
    float new_dE1=0.0;    
    float old_CC=0.0;
    float new_CC=0.0;
    float old_CC1=0.0;
    float new_CC1=0.0;    
    float old_dE2=0.0;
    float new_dE2=0.0; 
    float old_vwd1=0.0;
    float old_vwd2=0.0;
    float new_vwd1=0.0;
    float new_vwd2=0.0; 
    float new_dE=0.0;
    float old_dE=0.0;       
    float dE=0.0;
    float old_vwd= 0.0;
    float new_vwd= 0.0;
    float tmp_den=0.0;
    float tmp_rhc=0.0;
    float old_clash = 0.0;
    float new_clash = 0.0; 
    vector<float> coord_change(3,0.0);
    vector<poseCoord> best_model;
    vector<poseCoord > fin_maty;
    float best_E=10.0;
    vector<vector<float>> atm_idx(pnum,vector<float>(3,0.0));    

    ObjexxFCL::FArray3D< float > rhoC, rhoMask, rhoOmask, rhoO2mask;
    ObjexxFCL::FArray3D< std::complex<float> > FrhoC, FrhoMask, FrhoCmask, FrhoOmask, FrhoO2mask;
    ObjexxFCL::FArray3D< std::complex<float> > FrhoO, FrhoO2; 

    float effReso = std::max( 2.4+0.8*MRC_R , double(MRC_R ));
    float k=(M_PI/effReso)*(M_PI/effReso);
    float a=33.0;  // treat everything as ALA
    float C=a*pow(k/3.1415926,1.5);
    int tmp_i=0;
    ObjexxFCL::FArray3D< float > rhoC0;
    ObjexxFCL::FArray3D< float > inv_rho_mask0;
    rhoC0.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    inv_rho_mask0.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
        rhoC0[t]=0.0;
        inv_rho_mask0[t]=1.0;
    }
 
    string elt_i; 
    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;  

    // calc old_E
    bool old_TF = true;
    if(old_TF)
    {        
        vector<poseCoord>().swap(fin_maty);
        for(int j=0;j<pnum;j++)
        {
            fin_maty.push_back(pointsBx[3*j+0]);
            fin_maty.push_back(pointsBx[3*j+1]);
            fin_maty.push_back(pointsBx[3*j+2]);
        }
        // add side_chain center and CO atom
        vector<point3f> old_decstr_sg(pnum);
        for(int j=0;j<pnum;j++)
        {
            old_decstr_sg[j].ss2=numtoss(bb[j].sst);
            old_decstr_sg[j].ssm=numtoss(bb[j].sst);
            old_decstr_sg[j].stype=numtoss(bb[j].sst);
            old_decstr_sg[j].x=fin_maty[3*j+1].x_[0];
            old_decstr_sg[j].y=fin_maty[3*j+1].x_[1];
            old_decstr_sg[j].z=fin_maty[3*j+1].x_[2];
            old_decstr_sg[j].ptn.x=fin_maty[3*j+0].x_[0];
            old_decstr_sg[j].ptn.y=fin_maty[3*j+0].x_[1];
            old_decstr_sg[j].ptn.z=fin_maty[3*j+0].x_[2];
            old_decstr_sg[j].ptc.x=fin_maty[3*j+2].x_[0];
            old_decstr_sg[j].ptc.y=fin_maty[3*j+2].x_[1];
            old_decstr_sg[j].ptc.z=fin_maty[3*j+2].x_[2];
            string resng = pointsrenm[j];
            char resny = resng[0];
            old_decstr_sg[j].aaa = resny;   
            old_decstr_sg[j].iaa = (unsigned char) aminoid(resny);
            if(old_decstr_sg[j].iaa>19) old_decstr_sg[j].iaa=5;
            old_decstr_sg[j].resind=j+1;
//                tmp_tt=tmp_tt+1;
        }  
        str2tor(old_decstr_sg,pnum,3);
        tor2strsg2(old_decstr_sg,pnum,sgpos2);
        for(int j=0;j<pnum;j++)
        {
            // side_chain and CO
            fin_maty[3*j+1].x_sc[0]=old_decstr_sg[j].ptsg.x;
            fin_maty[3*j+1].x_sc[1]=old_decstr_sg[j].ptsg.y;
            fin_maty[3*j+1].x_sc[2]=old_decstr_sg[j].ptsg.z;
            fin_maty[3*j+2].x_o[0]=old_decstr_sg[j].pto.x;
            fin_maty[3*j+2].x_o[1]=old_decstr_sg[j].pto.y;
            fin_maty[3*j+2].x_o[2]=old_decstr_sg[j].pto.z;                        
        }  
        vector<float> del_ijx(3,0.0);
        vector<float> atm_jx(3,0.0);                
        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0); 
        tmp_i=0;
        for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
            rhoC0[t]=0.0;
            inv_rho_mask0[t]=1.0;
        }
    //        atm_idx=vector<vector<float> > ((rand_b2-rand_a2+1),vector<float>(3,0.0));
        atm_idx=vector<vector<float> > ((pnum),vector<float>(3,0.0));
        bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
        bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;
    //            for(int i=3*rand_a2;i<=3*rand_b2+2;i++)  // rotate fragment from beg_p to end_p
        for(int i=0;i<fin_maty.size();i++)  // rotate fragment from beg_p to end_p
        {
            if(fin_maty[i].elt_ == "CA")
            {           
                vector<float> cartX1;
                vector<float> fracX1;
                vector<float> atm_idg(3,0.0);
                elt_i = fin_maty[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMap.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_maty[i].x_;           
                MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idg[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                atm_idg[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                atm_idg[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMap.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idg[2]-atm_jx[2])/theDensityMap.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMap.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idg[1] - atm_jx[1])/theDensityMap.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMap.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idg[0] - atm_jx[0])/theDensityMap.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC0(x,y,z) += atm;
                            inv_rho_mask0(x,y,z) *= (1 - inv_msk);


                    //        if(x>bond0max) bond0max = x;
                    //        if(y>bond1max) bond1max = y;
                    //        if(z>bond2max) bond2max = z;
                    //        if(x<bond0min) bond0min = x;
                    //        if(y<bond1min) bond1min = y;
                    //        if(z<bond2min) bond2min = z;
    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                }
            } 
    /*        if(i%3==1)
            {
                vector<float> cartX1;
                vector<float> fracX1;
                vector<float> atm_idg(3,0.0);
                elt_i = fin_maty[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMap.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_maty[i].x_sc;           
                MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idg[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                atm_idg[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                atm_idg[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMap.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idg[2]-atm_jx[2])/theDensityMap.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMap.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idg[1] - atm_jx[1])/theDensityMap.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMap.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idg[0] - atm_jx[0])/theDensityMap.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC0(x,y,z) += atm;
                            inv_rho_mask0(x,y,z) *= (1 - inv_msk);


                    //        if(x>bond0max) bond0max = x;
                    //        if(y>bond1max) bond1max = y;
                    //        if(z>bond2max) bond2max = z;
                    //        if(x<bond0min) bond0min = x;
                    //        if(y<bond1min) bond1min = y;
                    //        if(z<bond2min) bond2min = z;
    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                }
            }
            if(i%3==2)
            {
                vector<float> cartX1;
                vector<float> fracX1;
                vector<float> atm_idg(3,0.0);
            //    elt_i = fin_maty[i].elt_;
                elt_i = 'O';
                elt_i = elt_i[0];
            //    cout<<elt_i<<endl;
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMap.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_maty[i].x_o;           
                MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idg[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                atm_idg[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                atm_idg[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMap.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idg[2]-atm_jx[2])/theDensityMap.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMap.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idg[1] - atm_jx[1])/theDensityMap.grid[1];
                        // wrap-around??
                        if ( del_ijx[1] > 0.5 ) del_ijx[1]-=1.0;
                        if ( del_ijx[1] < -0.5 ) del_ijx[1]+=1.0;                        
                        del_ijx[0] = 0.0;
                        vector<float> frac_tmpy;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpy);
                        if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                        
                        for(int x=1;x<=theDensityMap.density.u1();x++)
                        {
                            atm_jx[0] = x;
                            del_ijx[0] = (atm_idg[0] - atm_jx[0])/theDensityMap.grid[0];
                            // wrap-around??
                            if ( del_ijx[0] > 0.5 ) del_ijx[0]-=1.0;
                            if ( del_ijx[0] < -0.5 ) del_ijx[0]+=1.0;                            
                            vector<float> cart_del_ij2;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ijx,cart_del_ij2);
                            float d2 = square_len(cart_del_ij2);
                            if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                        
                            float atm = C*exp(-k*d2);
    //                            cout<<"MASK: "<<theDensityMap.ATOM_MASK<<endl;
                            float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK)  );
                        //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK)  );
                            float inv_msk = 1/(1+sigmoid_msk);
                            rhoC0(x,y,z) += atm;
                            inv_rho_mask0(x,y,z) *= (1 - inv_msk);


                //            if(x>bond0max) bond0max = x;
                //            if(y>bond1max) bond1max = y;
                //            if(z>bond2max) bond2max = z;
                //            if(x<bond0min) bond0min = x;
                //            if(y<bond1min) bond1min = y;
                //            if(z<bond2min) bond2min = z;
    //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
    //                                mask2(x,y,z) = 1.0; // problem?
    //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
    //                                    float atm = C*exp(-k*d2);
    //                                    rhoC(x,y,z) += atm;
    //                                }
    //                            }
                        }
                    }
                }
            } */
    //            cout<<" "<<tmp_i<<" "<<endl;           
            tmp_i = tmp_i + 1;
    //        } 
        } 
        float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
        float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
        float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
    //            tmp_den=0.0;
    //            tmp_rhc =0.0;
    //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
    //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
    //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
            // fetch this point
            clc_x2 = rhoC0[x];
    //            clc_x2 = tmp_rhc;
            obs_x2 = theDensityMap.density[x];
    //            obs_x2 = tmp_den;
            eps_x2 = 1-inv_rho_mask0[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
    //            eps_x2 = 1.0;

            // SMOOTHED
            sumCO_i2 += eps_x2*clc_x2*obs_x2;
            sumO_i2  += eps_x2*obs_x2;
            sumO2_i2 += eps_x2*obs_x2*obs_x2;
            sumC_i2  += eps_x2*clc_x2;
            sumC2_i2 += eps_x2*clc_x2*clc_x2;
            vol_i2   += eps_x2;
        }              
    //    mass_pdb[0] = mass_pdb[0]/allpdbmass;
    //    mass_pdb[1] = mass_pdb[1]/allpdbmass;
    //    mass_pdb[2] = mass_pdb[2]/allpdbmass;
    //    mass_mrc[0] = mass_mrc[0]/allmrcmass;
    //    mass_mrc[1] = mass_mrc[1]/allmrcmass;
    //    mass_mrc[2] = mass_mrc[2]/allmrcmass;
    //    olddis_p_m=Distance_point(mass_pdb,mass_mrc);     
        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
        if ( varC_i2 == 0 || varO_i2 == 0 ) {
            CC_i2 = 0;
        } else {
            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
        }   
        old_CC1=1.0-CC_i2;

//        vector<point3f>().swap(old_decstr_sg);
        float dis_p=0.0;
        old_clash = 0.0;
        for(int jj=0;jj<fin_maty.size();jj++)
        {
            if(fin_maty[jj].elt_=="CA")
            {
                string VR1 = fin_maty[jj].elt_;
                for(int t=jj;t<pointsBx.size();t++)
                {
                    if(pointsBx[t].elt_=="CA")
                    {
                    //    if(t>=3*rand_a && t<=3*rand_b+2) continue;
                        if(abs(t-jj)<3) continue;
                    //    if(abs(t)<3*end_p) continue;
                        string VR2 = fin_maty[t].elt_ ;
                        dis_p = Distance_point(fin_maty[jj].x_,fin_maty[t].x_);
                        old_clash = old_clash + GetVdwEgC(VR1,VR2,dis_p);

                //        cout<<"old_clash: "<<old_clash<<endl;
                //        if(dis_p < 3.70)
                //        {
                //            old_clash = old_clash + 1.0/ sqrt(dis_p + 0.0001) ;
                //            old_clash = old_clash + clash_s;
                        //    cout<<"jj t: "<<jj<<" "<<t<<" "<<old_clash<<endl;
                        //    cout<<"coor0: "<<fin_maty[jj].x_[0]<<" "<<fin_maty[jj].x_[1]<<" "<<fin_maty[jj].x_[2]<<endl;
                        //    cout<<"coor1: "<<pointsBx[t].x_[0]<<" "<<pointsBx[t].x_[1]<<" "<<pointsBx[t].x_[2]<<endl;
        //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                //        }
                    }
                }
            }
        }                    
        float old_Eres=0.0;
        for(int j=0;j<fin_maty.size();j++)
        {
            for(int js=j-3;js<=j+3;js++)
            {
                if(js==j) continue;
                if(js<0) continue;
                if(js>(pnum-1)) continue;
                float dx0 = Distance_point(fin_maty[j].x_,pointsBx0[js].x_);
        //        float dx1 = Distance_point(pointsBx0[j].x_,pointsBx0[js].x_);
        //        old_Eres = old_Eres + (dx0-dx1)*(dx0-dx1);
                old_Eres = old_Eres + dx0;
            }
        }                        
        old_Eres = old_Eres/float(pnum); 

        // Hbond energy
        float old_Ehbond=0.0;
        float hbalpha=0.0;
        float hbbeta=0.0;
    //    old_Ehbond=energyhbondnhoc2(old_decstr_sg,pnum,hbalpha,hbbeta);
        old_Ehbond=energyhbondnhoc2(old_decstr_sg,pnum);
    //    old_Ehbond = 5.0*old_Ehbond - hbalpha - hbbeta;

/*        // bondangle energy
        float old_fang=0.0;
        float old_Ebondang=0.0;
    //    cout<<"jjk: "<<endl;
        old_Ebondang =(float) Eenergybondangle(old_decstr_sg,pnum,old_fang);
    //    cout<<"jjk1: "<<endl;
        old_Ebondang = 0.30*old_Ebondang + old_fang; 
        old_Ebondang = old_Ebondang /float(pnum);        
        // bondlength energy                      
        float old_fene=0.0;
        float old_Ebondlen=0.0;
        old_Ebondlen = (float) energybondlength(old_decstr_sg,pnum,old_fene);
        old_Ebondlen = 0.5*old_Ebondlen + 20.0*old_fene;     

        // torsion angle
        float old_tor_fene=0.0;
        float old_tor_E=0.0;
        old_tor_E = energyrama(old_decstr_sg,pnum,old_tor_fene,rlogduke,ramaduke); 
        old_tor_E = 4.00*old_tor_E + old_tor_fene; */
        old_dE = 500.0*old_CC1 + old_clash+ 1.0*old_Eres + 2.0*old_Ehbond;//+ 1.0*old_tor_E + 0.5*old_Ebondlen;// + 0.2*old_tor_E ;//+ 10.0*old_Eres + 0.2*old_Ebondlen;// + old_Ebondang;// + 10.0*old_Eres; //// + old_Ebondang/1000.0;// + 200.0*(old_Ehbond/abs(pnum));
        cout<<"old_dE,old_CC1: "<<old_dE<<" "<<old_CC1<<endl;
    }


    cout<<"CCCCCC"<<endl;
    best_model = pointsBx;
    best_E = 100000.0;
    int NN0=0;
    int NN1=0;
    int NN2=0;
    int NN3=0;
    int NN4=0;
    int NN5=0;
    int NN6=0;
    int NN7=0;
    bool NN0TF=false;
    bool NN1TF=false;
    bool NN2TF=false;
    bool NN3TF=false;
    bool NN4TF=false;
    bool NN5TF=false;
    bool NN6TF=false;
    bool NN7TF=false;    
    // second structure rigid body movement
    bool calres = true;  
    vector<float> res_CC(pnum,0.0);
    vector<sssegment> segm;
    int iname =0;
    int movie=1;        
    int ssss=0;
    int algN=0;
    int rec_c = 5; // 5
    int rec_b = 10;
    int rec_a = 500;
    int rec_cb = rec_c*rec_b;
    int abc=0;       
    for(int jjy=0;jjy<rec_c;jjy++)
    {
        float KT0s = 1.0 ;//1.0;
        float KT0e = 0.5 ;//0.1;
        float KT0x= pow(float(KT0e/KT0s),float(float(jjy)/float(rec_c)));
        float KT0xx = KT0s*float(KT0x);
        float KTns = 0.1;
        float KTne = 0.05;
        float KTnx= pow(float(KTne/KTns),float(float(jjy)/float(rec_c)));
        float KTnxx = KTns*float(KTnx);  
        vector<float> Acprate(rec_b,0.0);
        for(int jjz=0;jjz<rec_b;jjz++)
        {
            int Aprate0=0;
            int Aprate1=0;            
//            float KTx= pow(float(KTnxx/KT0xx),float(float(jjz)/float(rec_b)));
//            float KT = KT0xx*float(KTx);  
            float KT = KTne + KT0s*(1.0-float(float(abc)/float(rec_cb)))*cos((1.0+float(abc)/float(float(rec_cb)/10.0))*PI*(float(abc)/float(rec_cb)))*cos((1.0+float(abc)/float(float(rec_cb)/10.0))*PI*(float(abc)/float(rec_cb)));
            abc = abc +1;            
            cout<<"KT: "<<KT<<endl; 
            vector<sssegment>().swap(segm);
            if(calres)
            {
                res_CC = vector<float>(pnum,0.0);
                for(int j=0;j<pnum;j++)
                {
                    vector<poseCoord> sing_res(1);
                    sing_res[0]= pointsBx[3*j+1];
                    res_CC[j] = theDensityMap.matchposex(sing_res);
            //        cout<<"j,CC: "<<j<<" "<<res_CC[j]<<endl; 
                }
                vector<float> res_CC3(pnum,0.0);
                for(int j=0;j<pnum;j++)
                {
                    for(int jp=-1;jp<2;jp++)
                    {
                        int jjp = j+ jp;
                        if(jjp<0) jjp=0;
                        if(jjp>pnum-1) jjp=pnum-1;
                        res_CC3[j] = res_CC3[j] + res_CC[jjp];
                    }
                    res_CC3[j] = res_CC3[j]/3.0;
                }
                for(int j=0;j<pnum;j++)
                {
                    int initx0=0;
                    int initx1=0;
                    if(res_CC3[j]<=0.1) // 0 is the cut off value for structure in the out of density map
                    {
                        int kx0=j-1;
                        if(kx0<0) kx0=0;
                        initx0=kx0;
                        for(int i=j;i<pnum;i++)
                        {
                            if(res_CC3[i]>0.1) 
                            {
                                float randcx=(2.0*randf0and1()-1.0);                                
                                if(randcx<0.0) break;
                        //        break;
                            }
                            j++;
                        }
                        int kx1=j+1;
                        if(kx1>pnum-1) kx1 = pnum-1;
                        initx1=kx1;
                    }
                    int dis_01 = initx1-initx0+1;
                    if(dis_01>=5) // 7 is the number which judgement how many residue should be as fragment
                    {
                        sssegment smg;
                        smg.init = initx0;
                        smg.term = initx1;
                        segm.push_back(smg);
                    }
                }
            }                   
            for(int jjx=0;jjx<rec_a;jjx++)
            {
                vector<poseCoord>().swap(fin_maty);
                for(int j=0;j<pnum;j++)
                {
                    fin_maty.push_back(pointsBx[3*j+0]);
                    fin_maty.push_back(pointsBx[3*j+1]);
                    fin_maty.push_back(pointsBx[3*j+2]);
                }              
                vector<sssegment> segmx;
                bool TF_segmx = false;
                if(segm.size()>0)
                {
                    TF_segmx=true;
                }
                if(TF_segmx)
                {
                    int i=0;
                    int last_segmx = segm.size()-1;
                    int ia1 = segm[i].init;
                    int ia2 = segm[i].term;
                    if(ia1<3)
                    {
                        sssegment smga;   
                        smga.init = ia1;
                        smga.term = ia2; 
                        for(int j=i+1;j<=last_segmx;j++)
                        {
                            int ic1 = segm[j].init;
                            int ic2 = segm[j].term;
                            if(abs(ic1-ia2)<7)
                            {
                                smga.term = ic2;
                                ia2 = ic2;
                                i++;
                            } else
                            {                
                                break;
                            }               
                        }                  
                        segmx.push_back(smga);                         
                    } else
                    {
                        sssegment smga;   
                        smga.init = ia1;
                        smga.term = ia2;  
                        segmx.push_back(smga);                       
                    }

                    int ib1 = segm[last_segmx].init;
                    int ib2 = segm[last_segmx].term; 
                    sssegment smgb;
                    if(ib2==(pnum-1))
                    {
            //            sssegment smg;   
                        smgb.init = ib1;
                        smgb.term = ib2; 
                        for(int j=last_segmx-1;j>=0;j--)
                        {
                            int ic1 = segm[j].init;
                            int ic2 = segm[j].term;
                            if(abs(ic2-ib1)<7)
                            {
                                smgb.init = ic1;
                                ib1 = ic1;
                                last_segmx--;
                    //            i++;
                            } else
                            {                
                                break;
                            }               
                        }                  
            //            segmx.push_back(smg);                                               
                    } else
                    {
                        smgb.init = ib1;
                        smgb.term = ib2;                        
                    }
                    for(int k=i+1;k<last_segmx;k++)
                    {
                        sssegment smg;
                        smg.init = segm[k].init;
                        smg.term = segm[k].term; 
                        segmx.push_back(smg);
                    }
                    segmx.push_back(smgb);
                }
                bool TF_segmy = false;
                if(segmx.size()>0)
                {
                    TF_segmy = true;
                }                   
                for(int jjk=0;jjk<1;jjk++)
                {

                    ang=90.0;
                    int rand_a=0,rand_b=0;
                    int rand_a2=0,rand_b2=0;
                    int randtx = randIntCustom(0,pnum-1);
                    { 
                        int randt = randIntCustom(0,4)+3;
                        int jmx=randtx,jnx=randtx;
                        float randt_num = 2.0*randf0and1()-1.0;
                        if(randt_num>0)
                        {
                            for(int j=0;j<randt;j++)
                            {
                                jnx=randtx+j;
                                if(jnx>=(pnum-1)) break;
                            }
                            rand_a=jmx;
                            rand_b=jnx;                                                
                        }else
                        {
                            for(int j=0;j<randt;j++)
                            {
                                jmx=randtx-j;
                                if(jmx<=0) break;
                            } 
                            rand_a=jmx;
                            rand_b=jnx;                                               
                        }                
                    }

                    if(TF_segmy)
                    {
                 //       float randt_vs = 2.0*randf0and1()-1.0;
                    //    if(randt_vs>=0)
                    //    {
                            for(int i=0;i<segmx.size();i++)
                            {
                                int io1 = segmx[i].init;
                                int io2 = segmx[i].term;
                             //   if(rand_b<io2 && rand_a>io1) 
                             //   {
                             //       rand_a = io1;
                             //       rand_b = io2;
                             //   }                           
                                if(rand_a<io2 && rand_a>io1) rand_a = io1;
                                if(rand_b<io2 && rand_b>io1) rand_b = io2;
                            }
                    //    } 
                    }
                    int NSS=0;
                    float randt_vg = 2.0*randf0and1()-1.0;
                    if(randt_vg>=0)
                    {
                        for(int i=1;i<pnum;i++)
                        {
                            int rand_ai = rand_a-i;
                            if(rand_ai<0) break;
                            if(numtoss(bb[rand_ai].sst)=='H'||numtoss(bb[rand_ai].sst)=='E')
                            {
                                rand_a = rand_ai;
                            } else
                            {
                                break;
                            }
                        }
                        for(int i=1;i<pnum;i++)
                        {
                            int rand_bi = rand_b+i;
                            if(rand_bi>(pnum-1)) break;
                        //    float randt_num = 2.0*randf0and1()-1.0;
                            if(numtoss(bb[rand_bi].sst)=='H'||numtoss(bb[rand_bi].sst)=='E')
                            {
                                rand_b = rand_bi;
                            } else
                            {
                                break;
                            }               
                        }                     
                    }  
                    rand_a2 = rand_a-2;
                    rand_b2 = rand_b+2;
                    if(rand_a <=0) rand_a = 0;
                    if(rand_b >(pnum-1)) rand_b= pnum-1;
                    if(rand_a2<=0) rand_a2=0;
                    if(rand_b2 >(pnum-1)) rand_b2= pnum-1;
                //    cout<<" "<<numtoss(bb[randtx].sst)<<endl;
                //    cout<<"rand a b: "<<rand_a<<" "<<rand_b<<endl;
                    if(abs(rand_b-rand_a+1)<2) continue;
                 //   cout<<"rand a b: "<<rand_a<<" "<<rand_b<<endl;
                    
                    
                    vector<float> coor_pdb(3,0.0);
                    coor_pdb=vector<float>(3,0.0);       

                    if(rand_b==(pnum-1)||rand_a==(0))
                    {
                        int rand_movement = randIntCustom(0,1000000)%5;
                        if(rand_movement == 0)
                        {
                            float asin_theta=2.0*randf0and1()-1.0;
                            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                            float apha=2.0*PI*randf0and1();
                            float awx=acos_theta*cos(apha);
                            float awy=acos_theta*sin(apha);
                            float awz=asin_theta;
                            // Translation Vector
                            float t0=0.0;
                            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                            float t3=(randf0and1()*2.0-1.0)*t0+0.0;
                        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

                            // Rotation matrix
                            float anggg=90.0;
                            float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                            float asin=sin(angle_rotategg*(PI/180.0));
                            float acos=cos(angle_rotategg*(PI/180.0));
                            float u[3][3];
                            u[0][0]=acos+awx*awx*(1.0-acos);
                            u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                            u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                            u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                            u[1][1]=acos+awy*awy*(1.0-acos);
                            u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                            u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                            u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                            u[2][2]=acos+awz*awz*(1.0-acos);
                            int rand_point = rand_a;
                            int tmp3=0;
                            vector<poseCoord > fin_mat;
                        //    fin_mat = tmp_mat;
                            tmp3=0;
                            vector<poseCoord>().swap(fin_mat);
                            fin_mat = fin_maty;
                        //    old_dE = old_clash + olddis_p_m;                      
                            vector<poseCoord > fin_matx=fin_mat;
                            tmp3=0;
                            if(rand_b==(pnum-1))
                            {
                                rand_point = rand_a;
                                axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                                axyz[1]=pointsBx[3*rand_point+1].x_[1];
                                axyz[2]=pointsBx[3*rand_point+1].x_[2];
                                for(int j=3*rand_a+1;j<=3*rand_b+2;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];                           
                                    tmp3=tmp3+1;
                                } 

                            } else
                            {
                                rand_point = rand_b;
                                axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                                axyz[1]=pointsBx[3*rand_point+1].x_[1];
                                axyz[2]=pointsBx[3*rand_point+1].x_[2];
                                for(int j=3*rand_a;j<=3*rand_b+1;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }
                                                         
                            }
                            vector<poseCoord>().swap(fin_maty);
                            fin_maty = fin_mat;                    
                            vector<poseCoord >().swap(fin_matx);
                        //    vector<poseCoord >().swap(fin_maty);
                        //    fin_maty = fin_mat;
                            vector<poseCoord >().swap(fin_mat);   
                            NN2TF = true ;                                                        
                        }
                        if(rand_movement == 1)
                        {
                            float asin_theta=2.0*randf0and1()-1.0;
                            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                            float apha=2.0*PI*randf0and1();
                            float awx=acos_theta*cos(apha);
                            float awy=acos_theta*sin(apha);
                            float awz=asin_theta;
                            // Translation Vector
                            float t0=0.3;
                            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                            float t3=(randf0and1()*2.0-1.0)*t0+0.0;
                        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

                            // Rotation matrix
                            float anggg=90.0;
                            float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                            float asin=sin(angle_rotategg*(PI/180.0));
                            float acos=cos(angle_rotategg*(PI/180.0));
                            float u[3][3];
                            u[0][0]=acos+awx*awx*(1.0-acos);
                            u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                            u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                            u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                            u[1][1]=acos+awy*awy*(1.0-acos);
                            u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                            u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                            u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                            u[2][2]=acos+awz*awz*(1.0-acos);            
                            // Rotation points
                    //      int rand_point = rand()%(pnumx);
                    //      int rand_point = rand()%(pnum);
                    //      axyz[0]=pointsx[rand_point].ca[0];
                    //      axyz[1]=pointsx[rand_point].ca[1];
                    //      axyz[2]=pointsx[rand_point].ca[2];
                        //    int rand_point = rand()%6+83 ;
                        //    int rand_point = randop ;
                        //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                        //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                        //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
                            coor_pdb=vector<float>(3,0.0);
                            for(int j=rand_a;j<=rand_b;j++)
                            {
                                coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                                coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                                coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];             
                            }
                            axyz[0]= coor_pdb[0]/float(rand_b-rand_a+1);  
                            axyz[1]= coor_pdb[1]/float(rand_b-rand_a+1);
                            axyz[2]= coor_pdb[2]/float(rand_b-rand_a+1);
                    //        cout<<"cent: "<< randtx<<" "<< axyz[0]<<" "<<axyz[1]<<" "<< axyz[2]<<endl; 

                        //    beg_p=rand_point;   

                            int tmp3=0;
                            vector<poseCoord > fin_mat;
                        //    fin_mat = tmp_mat;
                            tmp3=0;
                            vector<poseCoord>().swap(fin_mat);
                            fin_mat = fin_maty;
                        //    old_dE = old_clash + olddis_p_m;                      
                            vector<poseCoord > fin_matx=fin_mat;
                            tmp3=0;
                            for(int j=3*rand_a;j<=3*rand_b+2;j++) // not last N
                            {
                                fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                tmp3=tmp3+1;
                            }
                            vector<poseCoord >().swap(fin_matx);
                            int numseq=pnum;
                        //    if((rand_b-rand_a-3)<2) continue;
                        //    if(rand_a==0) continue;
                            vector<point3f> decstr(numseq);

                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=fin_mat[3*j+1].x_[0];
                                decstr[j].y=fin_mat[3*j+1].x_[1];
                                decstr[j].z=fin_mat[3*j+1].x_[2];
                                decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }
                            
                        //    if(rand_b+2>pnum-1) continue;
                            if(rand_b==(pnum-1))
                            {
                                if((rand_a2)<1||(rand_a+3)>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+3);
                            } else
                            {
                                if((rand_b-2)<1||rand_b2>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                            }
                            tmp_tt=0;
                            int tmp_ttx=0;
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*j+1].x_[0]=decstr[j].x;
                                fin_maty[3*j+1].x_[1]=decstr[j].y;
                                fin_maty[3*j+1].x_[2]=decstr[j].z;
                                fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                tmp_tt=tmp_tt+1;
                            }
                            vector<point3f>().swap(decstr);   
                            NN1TF = true;                                      
                        }
                        if(rand_movement==2)
                        {
                //            if((rand_b+1)>(pnum-1)) continue;
                            vector<float> sftv(3,0.0);
                            sftv[0] = (2.0*randf0and1()-1.0);
                            sftv[1] = (2.0*randf0and1()-1.0);
                            sftv[2] = (2.0*randf0and1()-1.0);
                        //    cout<<"sftv: "<<sftv[0]<<" "<<sftv[1]<<" "<<sftv[2]<<endl;
                            for(int j=rand_a;j<=rand_b;j++)
                            {
                                fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                            }
                            int numseq=pnum;
                        //    if((rand_b-rand_a-3)<2) continue;
                        //    if(rand_a==0) continue;
                            vector<point3f> decstr(numseq);

                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=fin_maty[3*j+1].x_[0];
                                decstr[j].y=fin_maty[3*j+1].x_[1];
                                decstr[j].z=fin_maty[3*j+1].x_[2];
                                decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }
                            if(rand_b==(pnum-1))
                            {
                                if((rand_a2)<1||(rand_a+3)>(pnum-1)) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+3);
                            } else
                            {
                                if((rand_b2)>(pnum-1)||(rand_b-2)<1) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2); 
                            }
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*j+1].x_[0]=decstr[j].x;
                                fin_maty[3*j+1].x_[1]=decstr[j].y;
                                fin_maty[3*j+1].x_[2]=decstr[j].z;
                                fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                            }     
                            vector<point3f>().swap(decstr);  
                            NN3TF = true;              
                        } 
                        if(rand_movement==3)
                        {
                //            if((rand_b+1)>(pnum-1)) continue;
                            float rand_mov0=(randf0and1()*2.0-1.0);
                            if(rand_mov0>=0.0)
                            {
                                vector<float> sftv(3,0.0);
                                sftv[0] = fin_maty[rand_b].x_[0] -fin_maty[rand_a].x_[0];
                                sftv[1] = fin_maty[rand_b].x_[1] -fin_maty[rand_a].x_[1];
                                sftv[2] = fin_maty[rand_b].x_[2] -fin_maty[rand_a].x_[2];
                                sftv[0] = sftv[0]/float(rand_b-rand_a+1);
                                sftv[1] = sftv[1]/float(rand_b-rand_a+1);
                                sftv[2] = sftv[2]/float(rand_b-rand_a+1);
                        //        cout<<"sftv: "<<sftv[0]<<" "<<sftv[1]<<" "<<sftv[2]<<endl;
    //                            sftv[0] = (2.0*randf0and1()-1.0);
    //                            sftv[1] = (2.0*randf0and1()-1.0);
    //                            sftv[2] = (2.0*randf0and1()-1.0);
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                    fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                    fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                    fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                    fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                    fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                    fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                    fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                    fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                                }
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=fin_maty[3*j+1].x_[0];
                                    decstr[j].y=fin_maty[3*j+1].x_[1];
                                    decstr[j].z=fin_maty[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                            //    if(rand_b>(pnum-3)) 
                            //    {
                            //        int rand_a_a = rand_a-2;
                            //        int rand_a_b = rand_a+2;
                            //        if(rand_a+2>pnum-1) rand_a_b = pnum-1;
                            //        if((rand_a-2)<1) rand_a_a = 1;
                            //        mcfragsweepLMP2(decstr,numseq,rand_a_a,rand_a_b);
                            //        mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+1);
                            //    } else
                            //    {
                            //        int rand_b_a = rand_b-2;
                            //        int rand_b_b = rand_b+2;
                            //        if((rand_b-2)<1) rand_b_a=1;
                            //        if(rand_b+2>pnum-1) rand_b_b = pnum-1;
                            //        mcfragsweepLMP2(decstr,numseq,rand_b_a,rand_b_b);
                            //        mcfragsweepLMP2(decstr,numseq,rand_b-1,rand_b2); 
                            //    }
                                if(rand_b==(pnum-1))
                                {
                                    if((rand_a2)<1||(rand_a+3)>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+3);
                                } else
                                {
                                    if((rand_b-2)<1||rand_b2>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                }                            
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*j+1].x_[0]=decstr[j].x;
                                    fin_maty[3*j+1].x_[1]=decstr[j].y;
                                    fin_maty[3*j+1].x_[2]=decstr[j].z;
                                    fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                }     
                                vector<point3f>().swap(decstr);
                            } else
                            {
                                vector<float> sftv(3,0.0);
                                sftv[0] = -fin_maty[rand_b].x_[0] + fin_maty[rand_a].x_[0];
                                sftv[1] = -fin_maty[rand_b].x_[1] + fin_maty[rand_a].x_[1];
                                sftv[2] = -fin_maty[rand_b].x_[2] + fin_maty[rand_a].x_[2];
                                sftv[0] = sftv[0]/float(rand_b-rand_a);
                                sftv[1] = sftv[1]/float(rand_b-rand_a);
                                sftv[2] = sftv[2]/float(rand_b-rand_a);
                        //        cout<<"sftv: "<<sftv[0]<<" "<<sftv[1]<<" "<<sftv[2]<<endl;
    //                            sftv[0] = (2.0*randf0and1()-1.0);
    //                            sftv[1] = (2.0*randf0and1()-1.0);
    //                            sftv[2] = (2.0*randf0and1()-1.0);
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                    fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                    fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                    fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                    fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                    fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                    fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                    fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                    fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                                }
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=fin_maty[3*j+1].x_[0];
                                    decstr[j].y=fin_maty[3*j+1].x_[1];
                                    decstr[j].z=fin_maty[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                            //    if(rand_b>(pnum-2)) 
                            //    {
                            //        int rand_a_a = rand_a-2;
                            //        int rand_a_b = rand_a+2;
                            //        if(rand_a+2>pnum-1) rand_a_b = pnum-1;
                            //        if((rand_a-2)<1) rand_a_a = 1;
                            //        mcfragsweepLMP2(decstr,numseq,rand_a_a,rand_a_b);
                            //        mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+1);
                            //    } else
                            //    {
                            //        int rand_b_a = rand_b-2;
                            //        int rand_b_b = rand_b+2;
                            //        if((rand_b-2)<1) rand_b_a=1;
                            //        if(rand_b+2>pnum-1) rand_b_b = pnum-1;
                            //        mcfragsweepLMP2(decstr,numseq,rand_b_a,rand_b_b);
                            //        mcfragsweepLMP2(decstr,numseq,rand_b-1,rand_b2); 
                            //    }
                                if(rand_b==(pnum-1))
                                {
                                    if((rand_a2)<1||(rand_a+3)>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+3);
                                } else
                                {
                                    if((rand_b-2)<1||rand_b2>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                }                            
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*j+1].x_[0]=decstr[j].x;
                                    fin_maty[3*j+1].x_[1]=decstr[j].y;
                                    fin_maty[3*j+1].x_[2]=decstr[j].z;
                                    fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                }     
                                vector<point3f>().swap(decstr);
                            }
                            NN4TF = true;
                        }
                        if(rand_movement==4)  // random change the coordinate from atom rand_a to atom atom rand_b
                        {
                //            int numseq=abs(rand_b-rand_a+1);
                            int numseq=pnum;
                    //        if((rand_b-rand_a-3)<2) continue;
                            if(rand_a==0) rand_a=1;
                            vector<point3f> decstr(numseq);

                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=pointsBx[3*j+1].x_[0];
                                decstr[j].y=pointsBx[3*j+1].x_[1];
                                decstr[j].z=pointsBx[3*j+1].x_[2];
                                decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                                decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                                decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                                decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                                decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                                decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }                               

                        //    mcfragsweepLMP2(decstr,numseq,rand_a,rand_b);
                            mcmovementLMP(decstr,numseq,rand_a,rand_b);
                            tmp_tt=0;
                            int tmp_ttx=0;
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                                fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                                fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                                fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                tmp_tt=tmp_tt+1;
                            }   
                            vector<point3f>().swap(decstr);  
                            NN6TF=true;           
                        }                                          
                    } else
                    {
                            int rand_movement = randIntCustom(0,10);
                       //     cout<<"sssssssssssssssss: "<<rand_movement<<endl;
                            if(rand_movement==0) // roatate around axis 
                            {
                                angle_rotate = (2.0*randf0and1()-1.0)*90.0; // rotate angle
                                GroupRotationpid(fin_maty[3*rand_a+1].x_,fin_maty[3*rand_b+2-1].x_,angle_rotate,fin_maty,3*rand_a+1,3*rand_b+2-1);  
                                NN0TF =true;                     
                            }
                            if(rand_movement==1)  // random change the coordinate from atom rand_a to atom atom rand_b
                            {
                    //            int numseq=abs(rand_b-rand_a+1);
                                int numseq=pnum;
                        //        if((rand_b-rand_a-3)<2) continue;
                        //        if(rand_a==0) rand_a=1;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=pointsBx[3*j+1].x_[0];
                                    decstr[j].y=pointsBx[3*j+1].x_[1];
                                    decstr[j].z=pointsBx[3*j+1].x_[2];
                                    decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                                    decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                                    decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                                    decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                                    decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                                    decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    decstr[j].aaa = resny;
                                    decstr[j].iaa = aminoid(resny);
                                    if(decstr[j].iaa>19) decstr[j].iaa=5;
                                    decstr[j].resind=j+1;                                    
                    //                tmp_tt=tmp_tt+1;
                                }                               
                                str2tor(decstr,pnum,3); 
                            //    mcfragsweepLMP2(decstr,numseq,rand_a,rand_b);
                                mcmovementLMP(decstr,numseq,rand_a,rand_b);
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                                    fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                                    fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                                    fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }   
                                vector<point3f>().swap(decstr);  
                                NN6TF=true;           
                            }
                            if(rand_movement == 2)
                            {     
                                int numseq=pnum;
                        //        if((rand_b-rand_a-3)<2) continue;
                        //        if(rand_a==0) rand_a=1;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=pointsBx[3*j+1].x_[0];
                                    decstr[j].y=pointsBx[3*j+1].x_[1];
                                    decstr[j].z=pointsBx[3*j+1].x_[2];
                                    decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                                    decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                                    decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                                    decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                                    decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                                    decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    decstr[j].aaa = resny;
                                    decstr[j].iaa = aminoid(resny);
                                    if(decstr[j].iaa>19) decstr[j].iaa=5;
                                    decstr[j].resind=j+1;                                    
                    //                tmp_tt=tmp_tt+1;
                                }                               
                                str2tor(decstr,pnum,3);                                                       
                                mcfragsweepome(decstr,pnum,phipsiprob);
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                                    fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                                    fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                                    fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }   
                                vector<point3f>().swap(decstr);                                
                            } 
                            if(rand_movement == 3)
                            {    
                                int numseq=pnum;
                        //        if((rand_b-rand_a-3)<2) continue;
                        //        if(rand_a==0) rand_a=1;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=pointsBx[3*j+1].x_[0];
                                    decstr[j].y=pointsBx[3*j+1].x_[1];
                                    decstr[j].z=pointsBx[3*j+1].x_[2];
                                    decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                                    decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                                    decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                                    decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                                    decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                                    decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    decstr[j].aaa = resny;
                                    decstr[j].iaa = aminoid(resny);
                                    if(decstr[j].iaa>19) decstr[j].iaa=5;
                                    decstr[j].resind=j+1;                                    
                    //                tmp_tt=tmp_tt+1;
                                }                               
                                str2tor(decstr,pnum,3);                                                       
                                mcfragsweepphi(decstr,pnum,phipsiprob);  
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                                    fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                                    fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                                    fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }   
                                vector<point3f>().swap(decstr);                                      
                            }     
                           if(rand_movement==4)
                            {      
                                int numseq=pnum;
                        //        if((rand_b-rand_a-3)<2) continue;
                        //        if(rand_a==0) rand_a=1;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=pointsBx[3*j+1].x_[0];
                                    decstr[j].y=pointsBx[3*j+1].x_[1];
                                    decstr[j].z=pointsBx[3*j+1].x_[2];
                                    decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                                    decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                                    decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                                    decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                                    decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                                    decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    decstr[j].aaa = resny;
                                    decstr[j].iaa = aminoid(resny);
                                    if(decstr[j].iaa>19) decstr[j].iaa=5;
                                    decstr[j].resind=j+1;                                    
                    //                tmp_tt=tmp_tt+1;
                                }                               
                                str2tor(decstr,pnum,3);                                                      
                                mcfragsweeppsi(decstr,pnum,phipsiprob);    
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                                    fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                                    fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                                    fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }   
                                vector<point3f>().swap(decstr);                                    
                            } 
                        if(rand_movement==5)
                        {
                //            bool flagok0=false;
                //            for(int j=rand_a;j<=rand_b;j++)
                //            {
                //                if(numtoss(bb[j].sst)!='C') flagok0 = true;
                //            }
                //            if(flagok0) continue; 
                //            if((rand_b+1)>(pnum-1)) continue;
                            vector<float> sftv(3,0.0);
                            sftv[0] = (2.0*randf0and1()-1.0);
                            sftv[1] = (2.0*randf0and1()-1.0);
                            sftv[2] = (2.0*randf0and1()-1.0);
                            for(int j=rand_a;j<=rand_b;j++)
                            {
                                fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                            }
                            int numseq=pnum;
                        //    if((rand_b-rand_a-3)<2) continue;
                        //    if(rand_a==0) continue;
                            vector<point3f> decstr(numseq);

                            for(int j=0;j<pnum;j++)
                            {
                                decstr[j].ss2=numtoss(bb[j].sst);
                                decstr[j].ssm=numtoss(bb[j].sst);
                                decstr[j].stype=numtoss(bb[j].sst);
                                decstr[j].x=fin_maty[3*j+1].x_[0];
                                decstr[j].y=fin_maty[3*j+1].x_[1];
                                decstr[j].z=fin_maty[3*j+1].x_[2];
                                decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                //                tmp_tt=tmp_tt+1;
                            }
                    //        if(randtx>(pnum-7))
                    //        {
                    //            if((rand_a-2)<1||(rand_a+2)>(pnum-1)) continue;
                                if(rand_a2<1) continue;
                                if(rand_a+2>(pnum-1)) continue;
                                if(rand_b-2<1) continue;
                                if(rand_b2>pnum-1) continue;                        
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                    //        } else
                    //        {
                    //            if((rand_b+2)>(pnum-1)||(rand_b-2)<1) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2); 
                    //        }
                            for(int j=0;j<pnum;j++)
                            {
                                fin_maty[3*j+1].x_[0]=decstr[j].x;
                                fin_maty[3*j+1].x_[1]=decstr[j].y;
                                fin_maty[3*j+1].x_[2]=decstr[j].z;
                                fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                            }     
                            vector<point3f>().swap(decstr);  
                            NN3TF = true;              
                        }
                            if(rand_movement == 6)
                            {
                        //        if(numtoss(bb[randtx].sst)=='E') continue;
                        //        cout<<"trdisx,anglex,KT: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;
                                float asin_theta=2.0*randf0and1()-1.0;
                                float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                                float apha=2.0*PI*randf0and1();
                                float awx=acos_theta*cos(apha);
                                float awy=acos_theta*sin(apha);
                                float awz=asin_theta;
                                // Translation Vector
                                float t0=0.3;
                                float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t3=(randf0and1()*2.0-1.0)*t0+0.0;
                            //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

                                // Rotation matrix
                                float anggg=90.0;
                                float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                                float asin=sin(angle_rotategg*(PI/180.0));
                                float acos=cos(angle_rotategg*(PI/180.0));
                                float u[3][3];
                                u[0][0]=acos+awx*awx*(1.0-acos);
                                u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                                u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                                u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                                u[1][1]=acos+awy*awy*(1.0-acos);
                                u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                                u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                                u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                                u[2][2]=acos+awz*awz*(1.0-acos);            
                                // Rotation points
                        //      int rand_point = rand()%(pnumx);
                        //      int rand_point = rand()%(pnum);
                        //      axyz[0]=pointsx[rand_point].ca[0];
                        //      axyz[1]=pointsx[rand_point].ca[1];
                        //      axyz[2]=pointsx[rand_point].ca[2];
                            //    int rand_point = rand()%6+83 ;
                            //    int rand_point = randop ;
                            //    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                            //    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                            //    axyz[2]=pointsBx[3*rand_point+1].x_[2];
                                coor_pdb=vector<float>(3,0.0);
                                for(int j=rand_a;j<=rand_b;j++)
                                {
                                    coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                                    coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                                    coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];
                                }                        
                                axyz[0]= coor_pdb[0]/float(rand_b-rand_a+1);  
                                axyz[1]= coor_pdb[1]/float(rand_b-rand_a+1);
                                axyz[2]= coor_pdb[2]/float(rand_b-rand_a+1);
                        //        cout<<"cent: "<< randtx<<" "<< axyz[0]<<" "<<axyz[1]<<" "<< axyz[2]<<endl; 

                            //    beg_p=rand_point;   

                                int tmp3=0;
                                vector<poseCoord > fin_mat;
                            //    fin_mat = tmp_mat;
                                tmp3=0;
                                vector<poseCoord>().swap(fin_mat);
                                fin_mat = fin_maty;
                    //            for(int j=rand_a+2;j<=rand_b-2;j++)
                    //            {
                    //                fin_mat.push_back(pointsBx[3*j+0]);
                    //                fin_mat.push_back(pointsBx[3*j+1]);
                    //                fin_mat.push_back(pointsBx[3*j+2]);               
                    //                tmp3=tmp3+1;
                    //            }        
                            //    old_dE = old_clash + olddis_p_m;                      
                                vector<poseCoord > fin_matx=fin_mat;
                                tmp3=0;
                                for(int j=3*rand_a;j<=3*rand_b+2;j++) // not last N
                                {
                                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                    tmp3=tmp3+1;
                                }
                                vector<poseCoord >().swap(fin_matx);
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=fin_mat[3*j+1].x_[0];
                                    decstr[j].y=fin_mat[3*j+1].x_[1];
                                    decstr[j].z=fin_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                                if(rand_a2<1) continue;
                                if(rand_a+2>(pnum-1)) continue;
                                if(rand_b-2<1) continue;
                                if(rand_b2>pnum-1) continue;
                                mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                                mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*j+1].x_[0]=decstr[j].x;
                                    fin_maty[3*j+1].x_[1]=decstr[j].y;
                                    fin_maty[3*j+1].x_[2]=decstr[j].z;
                                    fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }
                                vector<point3f>().swap(decstr); 
                                NN1TF=true;
                            } 
                            if(rand_movement == 7)
                            {
                                // Rotation axis 
                                float asin_theta=2.0*randf0and1()-1.0;
                                float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                                float apha=2.0*PI*randf0and1();
                                float awx=acos_theta*cos(apha);
                                float awy=acos_theta*sin(apha);
                                float awz=asin_theta;
                                // Translation Vector
                                float t0=0.0;
                                float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                                float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                                // Rotation matrix
                                float anggg=90.0;
                                float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                                float asin=sin(angle_rotategg*(PI/180.0));
                                float acos=cos(angle_rotategg*(PI/180.0));
                                float u[3][3];
                                u[0][0]=acos+awx*awx*(1.0-acos);
                                u[0][1]=awx*awy*(1.0-acos)-awz*asin;
                                u[0][2]=awx*awz*(1.0-acos)+awy*asin;
                                u[1][0]=awx*awy*(1.0-acos)+awz*asin;
                                u[1][1]=acos+awy*awy*(1.0-acos);
                                u[1][2]=awy*awz*(1.0-acos)-awx*asin;
                                u[2][0]=awx*awz*(1.0-acos)-awy*asin;
                                u[2][1]=awy*awz*(1.0-acos)+awx*asin;
                                u[2][2]=acos+awz*awz*(1.0-acos);            
                                // Rotation points
                        //      int rand_point = rand()%(pnumx);
                        //      int rand_point = rand()%(pnum);
                        //      axyz[0]=pointsx[rand_point].ca[0];
                        //      axyz[1]=pointsx[rand_point].ca[1];
                        //      axyz[2]=pointsx[rand_point].ca[2];
                        //        vector<float> coor_pdb(3,0.0);   
                        //        for(int j=rand_a;j<=rand_b;j++)
                        //        {
                        //            coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
                        //            coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
                        //            coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];             
                        //        }
                        //        coor_pdb[0] = coor_pdb[0]/abs(rand_b-rand_a+1);
                        //        coor_pdb[1] = coor_pdb[1]/abs(rand_b-rand_a+1);
                        //        coor_pdb[2] = coor_pdb[2]/abs(rand_b-rand_a+1);                                                   
                                int rand_point = rand_a;
                        //        axyz[0]= coor_pdb[0];  
                        //        axyz[1]= coor_pdb[1];
                        //        axyz[2]= coor_pdb[2];               
                                axyz[0] =pointsBx[3*rand_point+1].x_[0]; // CA atom
                                axyz[1] =pointsBx[3*rand_point+1].x_[1];
                                axyz[2] =pointsBx[3*rand_point+1].x_[2];
                                int tmp3=0;
                                vector<poseCoord > fin_mat;
                            //    fin_mat = tmp_mat;
                                vector<poseCoord>().swap(fin_mat);
                                fin_mat = fin_maty;
                          //      for(int j=rand_a+2;j<=rand_b-2;j++)
                            //    {
                            //        fin_mat.push_back(pointsBx[3*j+0]);
                            //        fin_mat.push_back(pointsBx[3*j+1]);
                            //        fin_mat.push_back(pointsBx[3*j+2]);               
                            //        tmp3=tmp3+1;
                            //    }        
                            //    old_dE = old_clash + olddis_p_m;                      
                                vector<poseCoord > fin_matx;
                                vector<poseCoord >().swap(fin_matx);
                                fin_matx=fin_mat;
                                tmp3=0;
                                float rand_mov0=(randf0and1()*2.0-1.0);
                                if(rand_mov0>0)
                                {
                                    rand_point = rand_a;
                                    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                                    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                                    axyz[2]=pointsBx[3*rand_point+1].x_[2];
                                    for(int j=3*rand_a+1;j<=3*rand_b+2;j++) // not last N
                                    {
                                        fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                        fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                        fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                        tmp3=tmp3+1;
                                    }                                                        
                                } else 
                                {
                                    rand_point = rand_b;
                                    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                                    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                                    axyz[2]=pointsBx[3*rand_point+1].x_[2];
                                    for(int j=3*rand_a;j<=3*rand_b+1;j++) // not last N
                                    {
                                        fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                                        fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                                        fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                                        tmp3=tmp3+1;
                                    }                            
                                }
                                vector<poseCoord >().swap(fin_matx);
                                int numseq=pnum;
                            //    if((rand_b-rand_a-3)<2) continue;
                            //    if(rand_a==0) continue;
                                vector<point3f> decstr(numseq);
                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=fin_mat[3*j+1].x_[0];
                                    decstr[j].y=fin_mat[3*j+1].x_[1];
                                    decstr[j].z=fin_mat[3*j+1].x_[2];
                                    decstr[j].ptn.x=fin_mat[3*j+0].x_[0];
                                    decstr[j].ptn.y=fin_mat[3*j+0].x_[1];
                                    decstr[j].ptn.z=fin_mat[3*j+0].x_[2];
                                    decstr[j].ptc.x=fin_mat[3*j+2].x_[0];
                                    decstr[j].ptc.y=fin_mat[3*j+2].x_[1];
                                    decstr[j].ptc.z=fin_mat[3*j+2].x_[2];
                    //                tmp_tt=tmp_tt+1;
                                }
                                if(rand_mov0>0)
                                {
                                    if((rand_b2)>(pnum-1)||(rand_b-1)<1) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_b-1,rand_b2);
                                } else
                                {
                                    if((rand_a2)<1||(rand_a+2)>(pnum-1)) continue;
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                                }
                                
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*j+1].x_[0]=decstr[j].x;
                                    fin_maty[3*j+1].x_[1]=decstr[j].y;
                                    fin_maty[3*j+1].x_[2]=decstr[j].z;
                                    fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }     
                                vector<point3f>().swap(decstr);  
                                NN2TF=true;              
                            }     
                           if(rand_movement==8)
                            {
                    //            bool flagok0=false;
                    //            for(int j=rand_a;j<=rand_b;j++)
                    //            {
                    //                if(numtoss(bb[j].sst)!='C') flagok0 = true;
                    //            }
                    //            if(flagok0) continue; 
                    //            if((rand_b+1)>(pnum-1)) continue;
                                float rand_mov0=(randf0and1()*2.0-1.0);
                                if(rand_mov0>=0.0)
                                {
                                    vector<float> sftv(3,0.0);
                                    sftv[0] = fin_maty[rand_b].x_[0] -fin_maty[rand_a].x_[0];
                                    sftv[1] = fin_maty[rand_b].x_[1] -fin_maty[rand_a].x_[1];
                                    sftv[2] = fin_maty[rand_b].x_[2] -fin_maty[rand_a].x_[2];
                                    sftv[0] = sftv[0]/float(rand_b-rand_a+1);
                                    sftv[1] = sftv[1]/float(rand_b-rand_a+1);
                                    sftv[2] = sftv[2]/float(rand_b-rand_a+1);
        //                            sftv[0] = (2.0*randf0and1()-1.0);
        //                            sftv[1] = (2.0*randf0and1()-1.0);
        //                            sftv[2] = (2.0*randf0and1()-1.0);
                                    for(int j=rand_a;j<=rand_b;j++)
                                    {
                                        fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                        fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                        fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                        fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                        fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                        fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                        fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                        fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                        fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                                    }
                                    int numseq=pnum;
                                //    if((rand_b-rand_a-3)<2) continue;
                                //    if(rand_a==0) continue;
                                    vector<point3f> decstr(numseq);

                                    for(int j=0;j<pnum;j++)
                                    {
                                        decstr[j].ss2=numtoss(bb[j].sst);
                                        decstr[j].ssm=numtoss(bb[j].sst);
                                        decstr[j].stype=numtoss(bb[j].sst);
                                        decstr[j].x=fin_maty[3*j+1].x_[0];
                                        decstr[j].y=fin_maty[3*j+1].x_[1];
                                        decstr[j].z=fin_maty[3*j+1].x_[2];
                                        decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                        decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                        decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                        decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                        decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                        decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                        //                tmp_tt=tmp_tt+1;
                                    }
                        ///            int rand_a_a = rand_a-2;
                        //            int rand_a_b = rand_a+2;
                        //            if(rand_a+2>pnum-1) rand_a_b = pnum-1;
                        //            if((rand_a-2)<1) rand_a_a = 1;
                        //            mcfragsweepLMP2(decstr,numseq,rand_a_a,rand_a_b);
                                    if(rand_a2<1) continue;
                                    if(rand_a+2>(pnum-1)) continue;
                                    if(rand_b-2<1) continue;
                                    if(rand_b2>pnum-1) continue;                                
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                        //            int rand_b_a = rand_b-2;
                        //            int rand_b_b = rand_b+2;
                        //            if((rand_b-2)<1) rand_b_a=1;
                        //            if(rand_b+2>pnum-1) rand_b_b = pnum-1;
                        //            mcfragsweepLMP2(decstr,numseq,rand_b_a,rand_b_b);
                                    mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                    for(int j=0;j<pnum;j++)
                                    {
                                        fin_maty[3*j+1].x_[0]=decstr[j].x;
                                        fin_maty[3*j+1].x_[1]=decstr[j].y;
                                        fin_maty[3*j+1].x_[2]=decstr[j].z;
                                        fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                        fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                        fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                        fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                        fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                        fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                    }     
                                    vector<point3f>().swap(decstr);
                                } else
                                {
                                    vector<float> sftv(3,0.0);
                                    sftv[0] = -fin_maty[rand_b].x_[0] + fin_maty[rand_a].x_[0];
                                    sftv[1] = -fin_maty[rand_b].x_[1] + fin_maty[rand_a].x_[1];
                                    sftv[2] = -fin_maty[rand_b].x_[2] + fin_maty[rand_a].x_[2];
                                    sftv[0] = sftv[0]/float(rand_b-rand_a);
                                    sftv[1] = sftv[1]/float(rand_b-rand_a);
                                    sftv[2] = sftv[2]/float(rand_b-rand_a);
        //                            sftv[0] = (2.0*randf0and1()-1.0);
        //                            sftv[1] = (2.0*randf0and1()-1.0);
        //                            sftv[2] = (2.0*randf0and1()-1.0);
                                    for(int j=rand_a;j<=rand_b;j++)
                                    {
                                        fin_maty[3*j+0].x_[0] = fin_maty[3*(j)+0].x_[0] + sftv[0];
                                        fin_maty[3*j+0].x_[1] = fin_maty[3*(j)+0].x_[1] + sftv[1];
                                        fin_maty[3*j+0].x_[2] = fin_maty[3*(j)+0].x_[2] + sftv[2];
                                        fin_maty[3*j+1].x_[0] = fin_maty[3*(j)+1].x_[0] + sftv[0];
                                        fin_maty[3*j+1].x_[1] = fin_maty[3*(j)+1].x_[1] + sftv[1];
                                        fin_maty[3*j+1].x_[2] = fin_maty[3*(j)+1].x_[2] + sftv[2];
                                        fin_maty[3*j+2].x_[0] = fin_maty[3*(j)+2].x_[0] + sftv[0];
                                        fin_maty[3*j+2].x_[1] = fin_maty[3*(j)+2].x_[1] + sftv[1];
                                        fin_maty[3*j+2].x_[2] = fin_maty[3*(j)+2].x_[2] + sftv[2];
                                    }
                                    int numseq=pnum;
                                //    if((rand_b-rand_a-3)<2) continue;
                                //    if(rand_a==0) continue;
                                    vector<point3f> decstr(numseq);

                                    for(int j=0;j<pnum;j++)
                                    {
                                        decstr[j].ss2=numtoss(bb[j].sst);
                                        decstr[j].ssm=numtoss(bb[j].sst);
                                        decstr[j].stype=numtoss(bb[j].sst);
                                        decstr[j].x=fin_maty[3*j+1].x_[0];
                                        decstr[j].y=fin_maty[3*j+1].x_[1];
                                        decstr[j].z=fin_maty[3*j+1].x_[2];
                                        decstr[j].ptn.x=fin_maty[3*j+0].x_[0];
                                        decstr[j].ptn.y=fin_maty[3*j+0].x_[1];
                                        decstr[j].ptn.z=fin_maty[3*j+0].x_[2];
                                        decstr[j].ptc.x=fin_maty[3*j+2].x_[0];
                                        decstr[j].ptc.y=fin_maty[3*j+2].x_[1];
                                        decstr[j].ptc.z=fin_maty[3*j+2].x_[2];
                        //                tmp_tt=tmp_tt+1;
                                    }
                            //        int rand_a_a = rand_a-2;
                            //        int rand_a_b = rand_a+2;
                            //        if(rand_a+2>pnum-1) rand_a_b = pnum-1;
                            //        if((rand_a-2)<1) rand_a_a = 1;
                            //        mcfragsweepLMP2(decstr,numseq,rand_a_a,rand_a_b);
                                    if(rand_a2<1) continue;
                                    if(rand_a+2>(pnum-1)) continue;
                                    if(rand_b-2<1) continue;
                                    if(rand_b2>pnum-1) continue;                                
                                    mcfragsweepLMP2(decstr,numseq,rand_a2,rand_a+2);
                            //        int rand_b_a = rand_b-2;
                            //        int rand_b_b = rand_b+2;
                            //        if((rand_b-2)<1) rand_b_a=1;
                            //        if(rand_b+2>pnum-1) rand_b_b = pnum-1;
                            //        mcfragsweepLMP2(decstr,numseq,rand_b_a,rand_b_b);
                                    mcfragsweepLMP2(decstr,numseq,rand_b-2,rand_b2);
                                    for(int j=0;j<pnum;j++)
                                    {
                                        fin_maty[3*j+1].x_[0]=decstr[j].x;
                                        fin_maty[3*j+1].x_[1]=decstr[j].y;
                                        fin_maty[3*j+1].x_[2]=decstr[j].z;
                                        fin_maty[3*j+0].x_[0]=decstr[j].ptn.x;
                                        fin_maty[3*j+0].x_[1]=decstr[j].ptn.y;
                                        fin_maty[3*j+0].x_[2]=decstr[j].ptn.z;
                                        fin_maty[3*j+2].x_[0]=decstr[j].ptc.x;
                                        fin_maty[3*j+2].x_[1]=decstr[j].ptc.y;
                                        fin_maty[3*j+2].x_[2]=decstr[j].ptc.z;
                                    }     
                                    vector<point3f>().swap(decstr);
                                }
                                NN4TF=true;
                            } 
                            if(rand_movement == 9)
                            {     
                                int numseq=pnum;
                        //        if((rand_b-rand_a-3)<2) continue;
                        //        if(rand_a==0) rand_a=1;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=pointsBx[3*j+1].x_[0];
                                    decstr[j].y=pointsBx[3*j+1].x_[1];
                                    decstr[j].z=pointsBx[3*j+1].x_[2];
                                    decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                                    decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                                    decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                                    decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                                    decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                                    decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    decstr[j].aaa = resny;
                                    decstr[j].iaa = aminoid(resny);
                                    if(decstr[j].iaa>19) decstr[j].iaa=5;
                                    decstr[j].resind=j+1;                                    
                    //                tmp_tt=tmp_tt+1;
                                }                               
                                str2tor(decstr,pnum,3);                                                       
                                mcfragsweeplen(decstr,pnum);
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                                    fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                                    fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                                    fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }   
                                vector<point3f>().swap(decstr);                                
                            } 
                            if(rand_movement == 10)
                            {     
                                int numseq=pnum;
                        //        if((rand_b-rand_a-3)<2) continue;
                        //        if(rand_a==0) rand_a=1;
                                vector<point3f> decstr(numseq);

                                for(int j=0;j<pnum;j++)
                                {
                                    decstr[j].ss2=numtoss(bb[j].sst);
                                    decstr[j].ssm=numtoss(bb[j].sst);
                                    decstr[j].stype=numtoss(bb[j].sst);
                                    decstr[j].x=pointsBx[3*j+1].x_[0];
                                    decstr[j].y=pointsBx[3*j+1].x_[1];
                                    decstr[j].z=pointsBx[3*j+1].x_[2];
                                    decstr[j].ptn.x=pointsBx[3*j+0].x_[0];
                                    decstr[j].ptn.y=pointsBx[3*j+0].x_[1];
                                    decstr[j].ptn.z=pointsBx[3*j+0].x_[2];
                                    decstr[j].ptc.x=pointsBx[3*j+2].x_[0];
                                    decstr[j].ptc.y=pointsBx[3*j+2].x_[1];
                                    decstr[j].ptc.z=pointsBx[3*j+2].x_[2];
                                    string resng = pointsrenm[j];
                                    char resny = resng[0];
                                    decstr[j].aaa = resny;
                                    decstr[j].iaa = aminoid(resny);
                                    if(decstr[j].iaa>19) decstr[j].iaa=5;
                                    decstr[j].resind=j+1;                                    
                    //                tmp_tt=tmp_tt+1;
                                }                               
                                str2tor(decstr,pnum,3);                                                       
                                mcfragsweepang(decstr,pnum);
                                tmp_tt=0;
                                int tmp_ttx=0;
                                for(int j=0;j<pnum;j++)
                                {
                                    fin_maty[3*tmp_tt+1].x_[0]=decstr[j].x;
                                    fin_maty[3*tmp_tt+1].x_[1]=decstr[j].y;
                                    fin_maty[3*tmp_tt+1].x_[2]=decstr[j].z;
                                    fin_maty[3*tmp_tt+0].x_[0]=decstr[j].ptn.x;
                                    fin_maty[3*tmp_tt+0].x_[1]=decstr[j].ptn.y;
                                    fin_maty[3*tmp_tt+0].x_[2]=decstr[j].ptn.z;
                                    fin_maty[3*tmp_tt+2].x_[0]=decstr[j].ptc.x;
                                    fin_maty[3*tmp_tt+2].x_[1]=decstr[j].ptc.y;
                                    fin_maty[3*tmp_tt+2].x_[2]=decstr[j].ptc.z;
                                    tmp_tt=tmp_tt+1;
                                }   
                                vector<point3f>().swap(decstr);                                
                            }                                                                                                     
                    }                    
                }   
        /*        if(algN>=50)
                {
                    mcsuperpose(pointsBx,theDensityMapx,pnum);
                    algN=0;
                }
                {
                    algN=algN+1;
                }  */
                

                vector<point3f> decstr_sg(pnum);
            //    vector<point3f>().swap(decstr_sg);
                for(int j=0;j<pnum;j++)
                {
                    decstr_sg[j].ss2=numtoss(bb[j].sst);
                    decstr_sg[j].ssm=numtoss(bb[j].sst);
                    decstr_sg[j].stype=numtoss(bb[j].sst);
                    decstr_sg[j].x=fin_maty[3*j+1].x_[0];
                    decstr_sg[j].y=fin_maty[3*j+1].x_[1];
                    decstr_sg[j].z=fin_maty[3*j+1].x_[2];
                    decstr_sg[j].ptn.x=fin_maty[3*j+0].x_[0];
                    decstr_sg[j].ptn.y=fin_maty[3*j+0].x_[1];
                    decstr_sg[j].ptn.z=fin_maty[3*j+0].x_[2];
                    decstr_sg[j].ptc.x=fin_maty[3*j+2].x_[0];
                    decstr_sg[j].ptc.y=fin_maty[3*j+2].x_[1];
                    decstr_sg[j].ptc.z=fin_maty[3*j+2].x_[2];
                    string resng = pointsrenm[j];
                    char resny = resng[0];
                    decstr_sg[j].aaa = resny;
                    decstr_sg[j].iaa = aminoid(resny);
                    if(decstr_sg[j].iaa>19) decstr_sg[j].iaa=5;
                    decstr_sg[j].resind=j+1;
    //                tmp_tt=tmp_tt+1;
                }  
                str2tor(decstr_sg,pnum,3);                  
        //        tor2strsg2(decstr_sg,pnum,sgpos2);
    /*            for(int j=0;j<pnum;j++)
                {
                    // side_chain and CO
                    fin_maty[3*j+1].x_sc[0]=decstr_sg[j].ptsg.x;
                    fin_maty[3*j+1].x_sc[1]=decstr_sg[j].ptsg.y;
                    fin_maty[3*j+1].x_sc[2]=decstr_sg[j].ptsg.z;
                    fin_maty[3*j+2].x_o[0]=decstr_sg[j].pto.x;
                    fin_maty[3*j+2].x_o[1]=decstr_sg[j].pto.y;
                    fin_maty[3*j+2].x_o[2]=decstr_sg[j].pto.z;                        
                }  */
                float new_Ehbond=0.0;
                float hbalpha=0.0;
                float hbbeta=0.0;
            //    new_Ehbond = energyhbondnhoc2(decstr_sg,pnum,hbalpha,hbbeta);
                new_Ehbond=energyhbondnhoc2(decstr_sg,pnum);
        //        new_Ehbond = 5.0*new_Ehbond - hbalpha - hbbeta;

                float new_Eres=0.0;
                for(int j=0;j<fin_maty.size();j++)
                {
                    for(int js=j-3;js<=j+3;js++)
                    {
                        if(js==j) continue;
                        if(js<0) continue;
                        if(js>(pnum-1)) continue;
                        float dx0 = Distance_point(fin_maty[j].x_,pointsBx0[js].x_);
        //                float dx1 = Distance_point(pointsBx0[j].x_,pointsBx0[js].x_);
        //                new_Eres = new_Eres + (dx0-dx1)*(dx0-dx1);
                        new_Eres = new_Eres + dx0;
                    }
                }                        
                new_Eres = new_Eres/float(pnum);   
                          
    /*            float new_fang=0.0;
                float new_Ebondang=0.0;
                new_Ebondang =(float) Eenergybondangle(decstr_sg,pnum,new_fang);
                new_Ebondang = 0.30*new_Ebondang + new_fang;  
                new_Ebondang = new_Ebondang/float(pnum); 

                float new_fene=0.0;
                float new_Ebondlen=0.0;
                new_Ebondlen = (float) energybondlength(decstr_sg,pnum,new_fene);
                new_Ebondlen = 0.5*new_Ebondlen + 20.0*new_fene;                  

                // torsion angle
                float new_tor_fene=0.0;
                float new_tor_E=0.0;
                new_tor_E = energyrama(decstr_sg,pnum,new_tor_fene,rlogduke,ramaduke); 
                new_tor_E = 4.00*new_tor_E + new_tor_fene;    */

                float dis_p=0.0;
                new_clash = 0.0;
                for(int jj=0;jj<fin_maty.size();jj++)
                {
                    if(fin_maty[jj].elt_=="CA")
                    {
            //            float VR1= GetVdwRadius(tmp_mat[jj].elt_);
                        string VR1 = fin_maty[jj].elt_;
                        for(int t=jj;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                //                if((jj+rand_a)==t) continue;
                        //        if(t>=(3*rand_a)&& t<=(3*rand_b+2)) continue;
                                if(abs(t-jj)<3) continue;
                //                if(t>3*rand_a && t<(3*rand_a+tmp_mat.size())) continue;
                //                float VR2= GetVdwRadius(pointsBx[t].elt_);
                                string VR2 = fin_maty[t].elt_ ;
                                dis_p = Distance_point(fin_maty[jj].x_,fin_maty[t].x_);
                                //if(dis_p < (VR1+VR2))
                        //        if(dis_p < 3.70)
                        //        {
                        //            new_clash = new_clash + 1.0/ sqrt(dis_p + 0.0001) ;
                 //                   old_clash = old_clash + exp( 3.75-dis_p) ;
                        //        }
                                  new_clash = new_clash + GetVdwEgC(VR1,VR2,dis_p);                        
                //                new_vwd = new_vwd + GetVdwEg(VR1,VR2,dis_p);
                //                new_vwd = new_vwd + GetVdwEgC(VR1,VR2,dis_p);
                //                if((dis_p)<2.0) tp_clashx = true;
                //                if(tp_clashx)
                //                {
                //                    cout<<"jj,t: "<<jj+rand_a<<" "<<t<<" dis: "<<VR1+VR2<<" "<<dis_p<<endl;
                //                    tp_clashx=false;
                //                }
            //                    cout<< "VR1+VR2,dis_p: "<<VR1+VR2<<" "<<dis_p<<endl;
                            }
                        }
        //                if(tp_clashx) break;
                    }
                }

                ObjexxFCL::FArray3D< float > rhoC2;
                ObjexxFCL::FArray3D< float > inv_rho_mask2;
                rhoC2.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
                inv_rho_mask2.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
                    for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
                        rhoC2[t]=0.0;
                        inv_rho_mask2[t]=1.0;
                    }
                vector<float> del_ij(3,0.0);
                vector<float> atm_j(3,0.0);             
                    del_ij=vector<float>(3,0.0);
                    atm_j=vector<float>(3,0.0);
                    int tm_i1 =0;
                    vector<vector<float>> atm_idy;
                    atm_idy=vector<vector<float> > ((pnum),vector<float>(3,0.0));
                bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
                bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;

//                    for(int j=3*rand_a2;j<=3*rand_b2+2;j++)
                for(int j=0;j<fin_maty.size();j++)
                {
                    if(fin_maty[j].elt_ == "CA")
                    {
                        vector<float> cartX1;
                        vector<float> fracX1;
                        vector<float> atm_idg(3,0.0); 
                    //    elt_i = tmp_mat[3*tm_i1+1].elt_;
                        elt_i = fin_maty[j].elt_;
                        elt_i = elt_i[0];
                        OneGaussianScattering sig_j = get_A( elt_i );
                        k = sig_j.k( theDensityMap.effectiveB );
                        C = sig_j.C( k );
                        if ( C < 1e-6 ) continue;

                //        cartX1 = tmp_mat[3*tm_i1+1].x_;
                        cartX1 = fin_maty[j].x_;
                //        tm_i1 = tm_i1 + 1;
                        MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                        atm_idg[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                        atm_idg[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                        atm_idg[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);                 
                        for(int z=1;z<=theDensityMap.density.u3();z++)
                        {
                            atm_j[2] = z;
                            del_ij[2] =(atm_idg[2]-atm_j[2])/theDensityMap.grid[2];
                            if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                            if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                            del_ij[0] = del_ij[1] = 0.0;
                            vector<float> frac_tmpz;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpz);
                            if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;               
                            for(int y=1;y<=theDensityMap.density.u2();y++)
                            {
                                atm_j[1] = y;
                                del_ij[1] = (atm_idg[1] - atm_j[1])/theDensityMap.grid[1];
                                // wrap-around??
                                if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                                if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                                del_ij[0] = 0.0;
                                vector<float> frac_tmpy;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpy);
                                if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                      
                                for(int x=1;x<=theDensityMap.density.u1();x++)
                                {
                                    atm_j[0] = x;
                                    del_ij[0] = (atm_idg[0] - atm_j[0])/theDensityMap.grid[0];
                                    // wrap-around??
                                    if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                    if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                    vector<float> cart_del_ij2;
                                    MatrixTimesTransVector(theDensityMap.f2c,del_ij,cart_del_ij2);
                                    float d2 = square_len(cart_del_ij2);
                                    if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                                
                                    float atm = C*exp(-k*d2);
                                    float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                    float inv_msk = 1.0/(1.0+sigmoid_msk);
                    //                rhoC2(x,y,z) += atm;
                    //                rhoC2(x,y,z) += atm;
                                    rhoC2(x,y,z) += atm;
                                    inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);


                //                    if(x>bond0max) bond0max = x;
                //                    if(y>bond1max) bond1max = y;
                //                    if(z>bond2max) bond2max = z;
                //                    if(x<bond0min) bond0min = x;
                //                    if(y<bond1min) bond1min = y;
                //                    if(z<bond2min) bond2min = z;
        //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
        //                                mask2(x,y,z) = 1.0; // problem?
        //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
        //                                    float atm = C*exp(-k*d2);
        //                                    rhoC(x,y,z) += atm;
        //                                }
        //                            }
                                }
                            }
                        }
                    }
        /*            if(j%3==1)
                    {
                        vector<float> cartX1;
                        vector<float> fracX1;
                        vector<float> atm_idg(3,0.0); 
                    //    elt_i = tmp_mat[3*tm_i1+1].elt_;
                        elt_i = fin_maty[j].elt_;
                        elt_i = elt_i[0];
                        OneGaussianScattering sig_j = get_A( elt_i );
                        k = sig_j.k( theDensityMap.effectiveB );
                        C = sig_j.C( k );
                        if ( C < 1e-6 ) continue;

                //        cartX1 = tmp_mat[3*tm_i1+1].x_;
                        cartX1 = fin_maty[j].x_sc;
                //        tm_i1 = tm_i1 + 1;
                        MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                        atm_idg[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                        atm_idg[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                        atm_idg[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);                 
                        for(int z=1;z<=theDensityMap.density.u3();z++)
                        {
                            atm_j[2] = z;
                            del_ij[2] =(atm_idg[2]-atm_j[2])/theDensityMap.grid[2];
                            if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                            if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                            del_ij[0] = del_ij[1] = 0.0;
                            vector<float> frac_tmpz;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpz);
                            if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;               
                            for(int y=1;y<=theDensityMap.density.u2();y++)
                            {
                                atm_j[1] = y;
                                del_ij[1] = (atm_idg[1] - atm_j[1])/theDensityMap.grid[1];
                                // wrap-around??
                                if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                                if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                                del_ij[0] = 0.0;
                                vector<float> frac_tmpy;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpy);
                                if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                      
                                for(int x=1;x<=theDensityMap.density.u1();x++)
                                {
                                    atm_j[0] = x;
                                    del_ij[0] = (atm_idg[0] - atm_j[0])/theDensityMap.grid[0];
                                    // wrap-around??
                                    if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                    if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                    vector<float> cart_del_ij2;
                                    MatrixTimesTransVector(theDensityMap.f2c,del_ij,cart_del_ij2);
                                    float d2 = square_len(cart_del_ij2);
                                    if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                                
                                    float atm = C*exp(-k*d2);
                                    float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                    float inv_msk = 1.0/(1.0+sigmoid_msk);
                    //                rhoC2(x,y,z) += atm;
                    //                rhoC2(x,y,z) += atm;
                                    rhoC2(x,y,z) += atm;
                                    inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);


            //                        if(x>bond0max) bond0max = x;
            //                        if(y>bond1max) bond1max = y;
            //                        if(z>bond2max) bond2max = z;
            //                        if(x<bond0min) bond0min = x;
            //                        if(y<bond1min) bond1min = y;
            //                        if(z<bond2min) bond2min = z;
        //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
        //                                mask2(x,y,z) = 1.0; // problem?
        //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
        //                                    float atm = C*exp(-k*d2);
        //                                    rhoC(x,y,z) += atm;
        //                                }
        //                            }
                                }
                            }
                        }
                    }
                    if(j%3==2)
                    {
                        vector<float> cartX1;
                        vector<float> fracX1;
                        vector<float> atm_idg(3,0.0); 
                    //    elt_i = tmp_mat[3*tm_i1+1].elt_;
                    //    elt_i = fin_maty[j].elt_;
                        elt_i = 'O';
                        elt_i = elt_i[0];
                        OneGaussianScattering sig_j = get_A( elt_i );
                        k = sig_j.k( theDensityMap.effectiveB );
                        C = sig_j.C( k );
                        if ( C < 1e-6 ) continue;

                //        cartX1 = tmp_mat[3*tm_i1+1].x_;
                        cartX1 = fin_maty[j].x_o;
                //        tm_i1 = tm_i1 + 1;
                        MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                        atm_idg[0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                        atm_idg[1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                        atm_idg[2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);                 
                        for(int z=1;z<=theDensityMap.density.u3();z++)
                        {
                            atm_j[2] = z;
                            del_ij[2] =(atm_idg[2]-atm_j[2])/theDensityMap.grid[2];
                            if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                            if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                            del_ij[0] = del_ij[1] = 0.0;
                            vector<float> frac_tmpz;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpz);
                            if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;               
                            for(int y=1;y<=theDensityMap.density.u2();y++)
                            {
                                atm_j[1] = y;
                                del_ij[1] = (atm_idg[1] - atm_j[1])/theDensityMap.grid[1];
                                // wrap-around??
                                if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                                if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                                del_ij[0] = 0.0;
                                vector<float> frac_tmpy;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpy);
                                if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                      
                                for(int x=1;x<=theDensityMap.density.u1();x++)
                                {
                                    atm_j[0] = x;
                                    del_ij[0] = (atm_idg[0] - atm_j[0])/theDensityMap.grid[0];
                                    // wrap-around??
                                    if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                    if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                    vector<float> cart_del_ij2;
                                    MatrixTimesTransVector(theDensityMap.f2c,del_ij,cart_del_ij2);
                                    float d2 = square_len(cart_del_ij2);
                                    if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;
                                
                                    float atm = C*exp(-k*d2);
                                    float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                                //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                    float inv_msk = 1.0/(1.0+sigmoid_msk);
                    //                rhoC2(x,y,z) += atm;
                    //                rhoC2(x,y,z) += atm;
                                    rhoC2(x,y,z) += atm;
                                    inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);


                //                    if(x>bond0max) bond0max = x;
                //                    if(y>bond1max) bond1max = y;
                //                    if(z>bond2max) bond2max = z;
                //                    if(x<bond0min) bond0min = x;
                //                    if(y<bond1min) bond1min = y;
                //                    if(z<bond2min) bond2min = z;
        //                            if ( d2 <= (theDensityMap.ATOM_MASK_SQ) ) {
        //                                mask2(x,y,z) = 1.0; // problem?
        //                                if ( d2 <= (theDensityMap.ATOM_DENS_SQ) ) {
        //                                    float atm = C*exp(-k*d2);
        //                                    rhoC(x,y,z) += atm;
        //                                }
        //                            }
                                }
                            }
                        }
                    }  */
                    tm_i1 = tm_i1 + 1;
                }
        //        rhoC2 = rhoC0;
        /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                    if(rhoC2[x]>max_rhc) max_rhc = rhoC2[x] ;
                }     */
        /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
            //        cout<<theDensityMap.density[x]<<" ";
                    if(rhoC2[x]<(3.0/5.0)*max_rhc) rhoC2[x]=0; 
                } */
    //            rhoC2 = rhoC0;  
                float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
                float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
                float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;                             
                sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
                sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
                clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
                for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
        //            tmp_den=0.0;
        //            tmp_rhc =0.0;
        //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
        //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
        //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                    // fetch this point
                    clc_x2 = rhoC2[x];
        //            clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMap.density[x];
        //            obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
        //            eps_x2 = 1.0;

                    // SMOOTHED
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                }
                varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
                varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
                if ( varC_i2 == 0 || varO_i2 == 0 ) {
                    CC_i2 = 0;
                } else {
                    CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
                } 
                new_CC1 =1.0-CC_i2;

                new_dE = 500.0*new_CC1 + new_clash + 1.0*new_Eres+ 2.0*new_Ehbond ;//+ 0.5*new_Ebondlen;//  + new_tor_E+ ;// + 0.2*new_tor_E;// + 10.0*new_Eres + 0.2*new_Ebondlen;//  + new_Ebondang ;//+ 10.0*new_Eres;//// + 1.0*Eres;// + new_Ebondang/1000.0;// + 200.0*(new_Ehbond/abs(pnum));
        //        
        //        cout<<"dE_CC, dE_clash,Hbond,bondlen,Eres,tor_E, bondang: "<< 500.0*(new_CC1-old_CC1)<<" "<<(new_clash-old_clash)<<" "<<1.0*(new_Ehbond-old_Ehbond)<<" "<<0.2*(new_Ebondlen-old_Ebondlen)<<" "<<10.0*(new_Eres-old_Eres)<<" "<<0.2*(new_tor_E-old_tor_E)<<" "<<(new_Ebondang-old_Ebondang )<<endl;
                cout<<"new_dE,new_CC1: "<<new_dE<<" "<<new_CC1<<endl;
                dE = new_dE - old_dE;
                if(new_dE <old_dE)
                {
                    old_dE = new_dE;
                    cout<<"KT,new_dE: "<<KT<<" "<<new_dE<<" "<<new_clash<<""<<new_Ehbond<<" "<<dE<<endl;
                //    cout<<"randx: "<<randtx<<endl;
                    int tmp_tmx=0;
                    for(int j=0;j<pnum;j++)
                    {
                        pointsBx[3*j+0] = fin_maty[3*j+0];
                        pointsBx[3*j+1] = fin_maty[3*j+1];
                        pointsBx[3*j+2] = fin_maty[3*j+2];
                        tmp_tmx = tmp_tmx +1;
                    }
                    if(new_dE<best_E)
                    {
                        best_E = new_dE;
                        vector<poseCoord>().swap(best_model);
                        best_model = fin_maty;
                    }  
                                vector<sssegment>().swap(segm);
                                if(calres)
                                {
                                    res_CC = vector<float>(pnum,2.0);
                                    for(int j=0;j<pnum;j++)
                                    {
                                        vector<poseCoord> sing_res(1);
                                        sing_res[0]= pointsBx[3*j+1];
                                        res_CC[j] = theDensityMap.matchposex(sing_res);
                                //        cout<<"j,CC: "<<j<<" "<<res_CC[j]<<endl; 
                                    }
                                    vector<float> res_CC3(pnum,0.0);
                                    for(int j=0;j<pnum;j++)
                                    {
                                        for(int jp=-1;jp<2;jp++)
                                        {
                                            int jjp = j+ jp;
                                            if(jjp<0) jjp=0;
                                            if(jjp>pnum-1) jjp=pnum-1;
                                            res_CC3[j] = res_CC3[j] + res_CC[jjp];
                                        }
                                        res_CC3[j] = res_CC3[j]/3.0;
                                    }                
                                    for(int j=0;j<pnum;j++)
                                    {
                                        int initx0=0;
                                        int initx1=0;
                                        if(res_CC3[j]<=0.05) // 0 is the cut off value for structure in the out of density map
                                        {
                                            int kx0=j-1;
                                            if(kx0<0) kx0=0;
                                            initx0=kx0;
                                            for(int i=j;i<pnum;i++)
                                            {
                                                if(res_CC3[i]>0.05) 
                                                {
                                                    float randcx=(2.0*randf0and1()-1.0);                                
                                                    if(randcx<0) break;
                                                }
                                                j++;
                                            }
                                            int kx1=j+1;
                                            if(kx1>pnum-1) kx1 = pnum-1;
                                            initx1=kx1;
                                        }
                                        int dis_01 = initx1-initx0+1;
                                        if(dis_01>=5) // 7 is the number which judgement how many residue should be as fragment
                                        {
                                            sssegment smg;
                                            smg.init = initx0;
                                            smg.term = initx1;
                                            segm.push_back(smg);
                                        }
                                    }
                                }                    
            /*            if(movie==1)
                        {
                            Model movie_model;
                            for(int i=0;i<posex.chains.size();i++)
                            {
                        //      cout<< "chainx size: "<<posex.chains.size()<<endl;
                        //      cout<< "chainy size: "<<posey.chains.size()<<endl;
                                Chain Chanx = posex.chains[i];
                                Chain Chany;
                                Chany.chaid = Chanx.chaid;
                            //    fin_pose.chains.push_back(Chanx);
                        //      Chain Chany = posey.chains[i];
                                int tmp_nx=0;
                                for(int j=0;j<pnum;j++) 
                                {
                        //          cout<<"j: "<<j<<endl;
                        //          cout<< "residuex size: "<< Chanx.residues.size()<<endl;
                        //          cout<< "residuey size: "<< Chany.residues.size()<<endl;
                                    Residue Resdx = Chanx.residues[j];
                                    Residue Resdy;
                                    Resdy.resseq = Resdx.resseq;
                                    Resdy.icode = Resdx.icode;
                                    Resdy.resname = Resdx.resname;
                                //    fin_pose.chains[i].residues.push_back(Resdx);
                        //          Residue Resdy = Chany.residues[j];
                        //          Rot rotx,roty;
                        //          cout<<"Res: "<<Resdx.atoms.size()<<endl;
                        //          cout<<"res: "<<points_pre3[tmp_nx].res_atm.size()<<endl;
                                    int tt=0;
                                    for(int t=0;t<Resdx.atoms.size();t++)
                                    {
                        //              cout<<"t: "<<t<<endl;
                        //              cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
                        //              cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
                                        Atom Atmx=Resdx.atoms[t];
                                        Atom Atmy;
                                    //    fin_pose.chains[i].residues[j].atoms.push_back(Atmx);
                        //              Atom Atmy=Resdy.heavyatoms[t];                                 
                                        if(Atmx.atona == " CA ") 
                                        {
                                //            fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[j*3+1].x_;
                                            Atmy = Atmx;
                                            Atmy.xyzVect = pointsBx[j*3+1].x_;
                                            Resdy.heavyatoms.push_back(Atmy);
                                            Resdy.atoms.push_back(Atmy);
                                    //        if(j==63) cout<<"xx: "<<pointsBx[j*3+1].x_[0]<<" "<<pointsBx[j*3+1].x_[1]<<" "<<pointsBx[j*3+1].x_[2]<<endl;
                                    //        fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[j].x_;

                                        } 
                                        if(Atmx.atona == " N  ") 
                                        {
                            //                fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[3*j].x_;
                                            Atmy = Atmx;
                                            Atmy.xyzVect = pointsBx[3*j].x_;
                                            Resdy.heavyatoms.push_back(Atmy);
                                            Resdy.atoms.push_back(Atmy);
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];
                                        }
                                        if(Atmx.atona == " C  ") 
                                        {
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[3*j+2].x_;
                                            Atmy = Atmx;
                                            Atmy.xyzVect = pointsBx[3*j+2].x_;
                                            Resdy.heavyatoms.push_back(Atmy);
                                            Resdy.atoms.push_back(Atmy);
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];
                                        }                                                    
                                    }   
                                    tmp_nx = tmp_nx + 1; 
                                    Chany.residues.push_back(Resdy);
                        //          cout<< "tmP: "<<tmp_nx<<endl;
                                }
                                movie_model.chains.push_back(Chany);
                            }
                            string movie_name;
                            movie_name = Int_to_String(iname)+ ".pdb";
                            char inputPDBxx[600];
                            writePDBStructure(inputPDBxx,movie_model,movie_name);
                            iname = iname +1;
                        }  */

            /*        for(int j=0;j<pnum;j++)
                    {
                        bb[j].indn = 3*j+0;
                        bb[j].indca = 3*j+1;
                        bb[j].indc = 3*j+2;  
                        string resng = pointsrenm[j];
                        char resny = resng[0];
                        bb[j].resid = aminoid(resny);
                        bb[j].sst = decstr_sg[j].ssm;                                                             
                    }     */
                //    calcsse2(bb, pnum, tmp_mat); 
                    calcssennhoc(bb,pnum,fin_maty);  
                    Aprate0 = Aprate0 + 1;  
                    if(NN0TF) NN0 = NN0+1;
                    if(NN1TF) NN1 = NN1+1;
                    if(NN2TF) NN2 = NN2+1;
                    if(NN3TF) NN3 = NN3+1;
                    if(NN4TF) NN4 = NN4+1;
                    if(NN5TF) NN5 = NN5+1;
                    if(NN6TF) NN6 = NN6+1;
                    if(NN7TF) NN7 = NN7+1;                      

                } else
                {
                    float tmpx=randf0and1();
                    float mc_v = exp(-dE/(KT)); 
                    if(tmpx < mc_v)
                    {
                        old_dE = new_dE;
                        int tmp_tmx=0;
                        for(int j=0;j<pnum;j++)
                        {
                            pointsBx[3*j+0] = fin_maty[3*j+0];
                            pointsBx[3*j+1] = fin_maty[3*j+1];
                            pointsBx[3*j+2] = fin_maty[3*j+2];
                        //    bb[j]=bby[tmp_tmx];
                            tmp_tmx = tmp_tmx +1;
                        }
                //        cout<<"KT,new_dE: "<<KT<<" "<<new_dE<<" "<<old_dE<<" "<<dE<<endl;
                        cout<<"KT,new_dE: "<<KT<<" "<<new_dE<<" "<<new_clash<<""<<new_Ehbond<<" "<<dE<<endl;

                 /*       for(int j=0;j<pnum;j++)
                        {
                            bb[j].indn = 3*j+0;
                            bb[j].indca = 3*j+1;
                            bb[j].indc = 3*j+2;  
                            string resng = pointsrenm[j];
                            char resny = resng[0];
                            bb[j].resid = aminoid(resny);
                            bb[j].sst = decstr_sg[j].ssm; 
                        }     */
                    //    calcsse2(bb, pnum, tmp_mat); 
                        calcssennhoc(bb,pnum,fin_maty); 
                                vector<sssegment>().swap(segm);
                                if(calres)
                                {
                                    res_CC = vector<float>(pnum,2.0);
                                    for(int j=0;j<pnum;j++)
                                    {
                                        vector<poseCoord> sing_res(1);
                                        sing_res[0]= pointsBx[3*j+1];
                                        res_CC[j] = theDensityMap.matchposex(sing_res);
                                //        cout<<"j,CC: "<<j<<" "<<res_CC[j]<<endl; 
                                    }
                                    vector<float> res_CC3(pnum,0.0);
                                    for(int j=0;j<pnum;j++)
                                    {
                                        for(int jp=-1;jp<2;jp++)
                                        {
                                            int jjp = j+ jp;
                                            if(jjp<0) jjp=0;
                                            if(jjp>pnum-1) jjp=pnum-1;
                                            res_CC3[j] = res_CC3[j] + res_CC[jjp];
                                        }
                                        res_CC3[j] = res_CC3[j]/3.0;
                                    }                
                                    for(int j=0;j<pnum;j++)
                                    {
                                        int initx0=0;
                                        int initx1=0;
                                        if(res_CC3[j]<=0.05) // 0 is the cut off value for structure in the out of density map
                                        {
                                            int kx0=j-1;
                                            if(kx0<0) kx0=0;
                                            initx0=kx0;
                                            for(int i=j;i<pnum;i++)
                                            {
                                                if(res_CC3[i]>0.05) 
                                                {
                                                    float randcx=(2.0*randf0and1()-1.0);                                
                                                    if(randcx<0) break;
                                                }
                                                j++;
                                            }
                                            int kx1=j+1;
                                            if(kx1>pnum-1) kx1 = pnum-1;
                                            initx1=kx1;
                                        }
                                        int dis_01 = initx1-initx0+1;
                                        if(dis_01>=5) // 7 is the number which judgement how many residue should be as fragment
                                        {
                                            sssegment smg;
                                            smg.init = initx0;
                                            smg.term = initx1;
                                            segm.push_back(smg);
                                        }
                                    }
                                }
                /*        if(movie==1)
                        {
                            Model movie_model;
                            for(int i=0;i<posex.chains.size();i++)
                            {
                        //      cout<< "chainx size: "<<posex.chains.size()<<endl;
                        //      cout<< "chainy size: "<<posey.chains.size()<<endl;
                                Chain Chanx = posex.chains[i];
                                Chain Chany;
                                Chany.chaid = Chanx.chaid;
                            //    fin_pose.chains.push_back(Chanx);
                        //      Chain Chany = posey.chains[i];
                                int tmp_nx=0;
                                for(int j=0;j<pnum;j++) 
                                {
                        //          cout<<"j: "<<j<<endl;
                        //          cout<< "residuex size: "<< Chanx.residues.size()<<endl;
                        //          cout<< "residuey size: "<< Chany.residues.size()<<endl;
                                    Residue Resdx = Chanx.residues[j];
                                    Residue Resdy;
                                    Resdy.resseq = Resdx.resseq;
                                    Resdy.icode = Resdx.icode;
                                    Resdy.resname = Resdx.resname;
                                //    fin_pose.chains[i].residues.push_back(Resdx);
                        //          Residue Resdy = Chany.residues[j];
                        //          Rot rotx,roty;
                        //          cout<<"Res: "<<Resdx.atoms.size()<<endl;
                        //          cout<<"res: "<<points_pre3[tmp_nx].res_atm.size()<<endl;
                                    int tt=0;
                                    for(int t=0;t<Resdx.atoms.size();t++)
                                    {
                        //              cout<<"t: "<<t<<endl;
                        //              cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
                        //              cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
                                        Atom Atmx=Resdx.atoms[t];
                                        Atom Atmy;
                                    //    fin_pose.chains[i].residues[j].atoms.push_back(Atmx);
                        //              Atom Atmy=Resdy.heavyatoms[t];                                 
                                        if(Atmx.atona == " CA ") 
                                        {
                                //            fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[j*3+1].x_;
                                            Atmy = Atmx;
                                            Atmy.xyzVect = pointsBx[j*3+1].x_;
                                            Resdy.heavyatoms.push_back(Atmy);
                                            Resdy.atoms.push_back(Atmy);
                                    //        if(j==63) cout<<"xx: "<<pointsBx[j*3+1].x_[0]<<" "<<pointsBx[j*3+1].x_[1]<<" "<<pointsBx[j*3+1].x_[2]<<endl;
                                    //        fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[j].x_;

                                        } 
                                        if(Atmx.atona == " N  ") 
                                        {
                            //                fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[3*j].x_;
                                            Atmy = Atmx;
                                            Atmy.xyzVect = pointsBx[3*j].x_;
                                            Resdy.heavyatoms.push_back(Atmy);
                                            Resdy.atoms.push_back(Atmy);
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];
                                        }
                                        if(Atmx.atona == " C  ") 
                                        {
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[3*j+2].x_;
                                            Atmy = Atmx;
                                            Atmy.xyzVect = pointsBx[3*j+2].x_;
                                            Resdy.heavyatoms.push_back(Atmy);
                                            Resdy.atoms.push_back(Atmy);
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
                        //                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];
                                        }                                                   
                                    }   
                                    tmp_nx = tmp_nx + 1; 
                                    Chany.residues.push_back(Resdy);
                        //          cout<< "tmP: "<<tmp_nx<<endl;
                                }
                                movie_model.chains.push_back(Chany);
                            }
                            string movie_name;
                            movie_name = Int_to_String(iname)+ ".pdb";
                            char inputPDBxx[600];
                            writePDBStructure(inputPDBxx,movie_model,movie_name);
                            iname = iname +1;
                        }  */
                        Aprate0 = Aprate0 + 1;
                        if(NN0TF) NN0 = NN0+1;
                        if(NN1TF) NN1 = NN1+1;
                        if(NN2TF) NN2 = NN2+1;
                        if(NN3TF) NN3 = NN3+1;
                        if(NN4TF) NN4 = NN4+1;
                        if(NN5TF) NN5 = NN5+1;
                        if(NN6TF) NN6 = NN6+1;
                        if(NN7TF) NN7 = NN7+1;                                               
                    }                   
                }
                NN0TF = false;
                NN1TF = false; 
                NN2TF = false; 
                NN3TF = false; 
                NN4TF = false; 
                NN5TF = false; 
                NN6TF = false;
                NN7TF = false;                

            }
            Acprate[jjz] = float(Aprate0)/float(rec_a);
        } 
        cout<<"ACP rate: ";
        for(int j=0;j<Acprate.size();j++)
        {
            cout<<Acprate[j]<<" ";
        //    if(j%15==0)
        //    {
        //        cout<<endl;
        //    }
        }
        cout<<endl;           
    } 
    for(int j=0;j<pnum;j++)
    {
        pointsBx[3*j+0] = best_model[3*j+0];
        pointsBx[3*j+1] = best_model[3*j+1];
        pointsBx[3*j+2] = best_model[3*j+2];
    //    bb[j]=bby[tmp_tmx];
    } 
    cout<<"NN0,NN1,NN2,NN3,NN4,NN5,NN6: "<<NN0<<" "<<NN1<<" "<<NN2<<" "<<NN3<<" "<<NN4<<" "<<NN5<<" "<<NN6<<endl;


    for(int i=0;i<posex.chains.size();i++)
    {
//      cout<< "chainx size: "<<posex.chains.size()<<endl;
//      cout<< "chainy size: "<<posey.chains.size()<<endl;
        Chain Chanx = posex.chains[i];
        Chain Chany;
        Chany.chaid = Chanx.chaid;
    //    fin_pose.chains.push_back(Chanx);
//      Chain Chany = posey.chains[i];
        int tmp_nx=0;
        for(int j=0;j<pnum;j++) 
        {
//          cout<<"j: "<<j<<endl;
//          cout<< "residuex size: "<< Chanx.residues.size()<<endl;
//          cout<< "residuey size: "<< Chany.residues.size()<<endl;
            Residue Resdx = Chanx.residues[j];
            Residue Resdy;
            Resdy.resseq = Resdx.resseq;
            Resdy.icode = Resdx.icode;
            Resdy.resname = Resdx.resname;
        //    fin_pose.chains[i].residues.push_back(Resdx);
//          Residue Resdy = Chany.residues[j];
//          Rot rotx,roty;
//          cout<<"Res: "<<Resdx.atoms.size()<<endl;
//          cout<<"res: "<<points_pre3[tmp_nx].res_atm.size()<<endl;
            int tt=0;
            for(int t=0;t<Resdx.atoms.size();t++)
            {
//              cout<<"t: "<<t<<endl;
//              cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
//              cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
                Atom Atmx=Resdx.atoms[t];
                Atom Atmy;
            //    fin_pose.chains[i].residues[j].atoms.push_back(Atmx);
//              Atom Atmy=Resdy.heavyatoms[t];                                 
                if(Atmx.atona == " CA ") 
                {
        //            fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[j*3+1].x_;
                    Atmy = Atmx;
                    Atmy.xyzVect = pointsBx[j*3+1].x_;
                    Resdy.heavyatoms.push_back(Atmy);
                    Resdy.atoms.push_back(Atmy);
            //        if(j==63) cout<<"xx: "<<pointsBx[j*3+1].x_[0]<<" "<<pointsBx[j*3+1].x_[1]<<" "<<pointsBx[j*3+1].x_[2]<<endl;
            //        fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[j].x_;

                } 
                if(Atmx.atona == " N  ") 
                {
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[3*j].x_;
                    Atmy = Atmx;
                    Atmy.xyzVect = pointsBx[3*j].x_;
                    Resdy.heavyatoms.push_back(Atmy);
                    Resdy.atoms.push_back(Atmy);
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];
                }
                if(Atmx.atona == " C  ") 
                {
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect = pointsBx[3*j+2].x_;
                    Atmy = Atmx;
                    Atmy.xyzVect = pointsBx[3*j+2].x_;
                    Resdy.heavyatoms.push_back(Atmy);
                    Resdy.atoms.push_back(Atmy);
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
//                    fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];
                }
        /*        if(Atmx.atona == " O  ") 
                {
    //              fin_pose.chains[i].residues[j].atoms[t].xyzVect = CBp[j].Ocd;
                    float C_O_length = 1.229;
                    float CA_C_O_angle = 120.09;
                    float N_CA_C_O_diangle = 179.67;
                    vector<float> O_coord;
                    O_coord=calculateCoordinates(pointsBx[3*j].x_, pointsBx[j*3+1].x_, pointsBx[3*j+2].x_, C_O_length, CA_C_O_angle, N_CA_C_O_diangle);
                    Atmy = Atmx;
                    Atmy.xyzVect = O_coord;
                    Resdy.heavyatoms.push_back(Atmy);
                    Resdy.atoms.push_back(Atmy);
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[3*j+2].x_[0] - pointsBx0[3*j+2].x_[0];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[3*j+2].x_[1] - pointsBx0[3*j+2].x_[1];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[3*j+2].x_[2] - pointsBx0[3*j+2].x_[2];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[0] + pointsBx[j].x_[0] - pointsBx0[j].x_[0];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[1] + pointsBx[j].x_[1] - pointsBx0[j].x_[1];
    //                fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] = fin_pose.chains[i].residues[j].atoms[t].xyzVect[2] + pointsBx[j].x_[2] - pointsBx0[j].x_[2];                
                }          */                                                     
            }   
            tmp_nx = tmp_nx + 1; 
            Chany.residues.push_back(Resdy);
//          cout<< "tmP: "<<tmp_nx<<endl;
        }
        fin_pose.chains.push_back(Chany);
    }

    cout<< "XXXXXXXXXX"<<endl;
    return fin_pose;            
}

int main(int argc, char* argv[])
{
    char inputPDB1[600];
    Model pose;
//  readPDB
    strcpy(inputPDB1,argv[1]);
    readpdbstructurex(inputPDB1,pose);
    cout<< "Read PDB file : "<<inputPDB1<<endl;

    string inputMRC;
    inputMRC = argv[2];

    float MRC_reso;
    float mapsampling = 0.0;
    MRC_reso = atof(argv[3]);
    mapsampling = atof(argv[4]);

    char bindir[600];
    strcpy(bindir,argv[5]); 

    clock_t start,finish;
    double totaltime;
    start=clock();         

    Model fin_p,fin_px;
    fin_p=REMC_sampley(pose,inputMRC,MRC_reso,mapsampling,bindir);

    string outnamex;
    vector<string> tmp_str,tmp_strx;
    tmp_str = string_splitx(argv[1],'.');
    tmp_strx = string_splitx(tmp_str[0],'/'); 
    string tmp_rx;
    if(tmp_strx.size()>0)
    {
        outnamex = tmp_strx[tmp_strx.size()-1] + "_11.pdb";
        tmp_rx = tmp_strx[tmp_strx.size()-1];
    }
    else
    {
        outnamex = tmp_strx[0] + "_11.pdb";
        tmp_rx = tmp_strx[0];
    }

    writePDBStructure(inputPDB1,fin_p,outnamex);
    char cmd[1000];
    sprintf(cmd, "%s/pulchra -c %s",bindir,outnamex.c_str());
    system(cmd);

//    sleep(10);

    char rmcmd[100];
    sprintf(rmcmd, "rm -f %s",outnamex.c_str());
    system(rmcmd);

    string inputPDB2=tmp_rx+ "_11.rebuilt.pdb";
    Model poseyy;
    readpdbstructurex(inputPDB2.data(),poseyy);   
    sprintf(rmcmd, "rm -f %s",inputPDB2.c_str());
    system(rmcmd);    

    string outname;
    if(tmp_strx.size()>0)
    {
        outname = tmp_strx[tmp_strx.size()-1] + "_outp00x03.pdb";
    }
    else
    {
        outname = tmp_strx[0] + "_outp00x03.pdb";
    }    

    cout<<"out: "<<outname<<endl;
    writePDBStructure(inputPDB1,poseyy,outname);

    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"\nrunning time: "<<totaltime<<" s！"<<endl;  
    return 0;    
}