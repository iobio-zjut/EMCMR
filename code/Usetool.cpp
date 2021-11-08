#include "MC.h"
#include "GetVdwRadius.h"
#include <time.h>
#include "Rotbuilder.h"
#include "CaDensityx.h"
#include "randomx.h"

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

Model REMC_sampley(Model posex, string inputM,float MRC_R,float mapsp,int dm_beg,int dm_end)
{
    int nresobins = 200;  //  # resolution bins for statistics
    vector<float> resobins, mapI, maskedmapI, modelI, maskI;
    vector<int> resobin_counts;
    vector<float> perResCC, perResStrain;

    vector<float> mapmapFSC, maskedMapMapFSC;
    vector<float> modelmapFSC, maskedModelMapFSC;

    bool truncate_map = false ;
    bool maskonly = false ;
    bool cutonly = false ;
    bool bin_squared =  false;
    bool mask =false ;
    string alt_mapfile=" ";

    // resolution limits for analysis
    float hires = 0.0;  // high res limit
    float lowres = 1000.0;  // low res limit
    float truncate_hires = 0.0;  // high res truncation
    float truncate_lowres = 1000.0;  // low res truncation
    float mask_resolution = 0.0;   //  radius for masking
    bool perres =false;  //dump extra output

    ObjexxFCL::FArray3D< float > rhoC, rhoMask, rhoOmask, rhoO2mask;
    ObjexxFCL::FArray3D< std::complex<float> > FrhoC, FrhoMask, FrhoCmask, FrhoOmask, FrhoO2mask;
    ObjexxFCL::FArray3D< std::complex<float> > FrhoO, FrhoO2;   

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
//  theDensityMap.calcRhoC(posey,4.0,rhoC,rhoMask);
//    theDensityMap.calcRhoCx(posex,4.0,rhoC,rhoMask);
//  theDensityMap.calcRhoCy(pose,4.0,rhoC,rhoMask);
//    float CC;
//  CC=theDensityMap.getRSCC(rhoC,rhoMask);
 //   CC=theDensityMap.getRSCCX(rhoC,rhoMask);
 //   cout<<"CC: "<<CC<<endl; 
    cout<<"origin: "<<theDensityMap.origin[0]<<" "<<theDensityMap.origin[1]<<" "<<theDensityMap.origin[2]<<endl;;

    Model fin_pose;
//    fin_pose = posex;
//  vector<Rot> fin_vec;
//  vector<Rot> beg_vec;
//  float min_dis=100.0;
//  vector<vector<float> > point;
//  vector<float> old_v(3,0.0),old_v1(3,0.0);
//  vector<vector<float> > old_v,old_v1,old_tmp,new_v;
//  old_v=vector<vector<float> >(4,vector<float>(3,0.0));
//  old_v1=vector<vector<float> >(4,vector<float>(3,0.0));
//  old_tmp = vector<vector<float> >(4,vector<float>(3,0.0));
//  new_v = vector<vector<float> >(4,vector<float>(3,0.0));
//  vector<float> new_v(3,0.0);
    float KT=0.001;
//  vector<vector<Rot> > frag_pre(7),frag_nat(7);


    
    vector<Rot> pointsx;
//  vector<vector<float> > pointsABx; // all backbone atom
    vector<poseCoord> pointsBx;
//  vector<vector<float> > pointsABx0; 
    vector<poseCoord> pointsBx0;
    vector<string> pointsrenm;
//  vector<vector<float> > pointsx,pointsy;
//  int pnum = posex.size();

//  cout<<"33 "<<endl;
//  cout<< posex.size()<<endl;
//  cout<< posey.size()<<endl;
//  for(int j=0;j<pnum;j++)
//  {
//      pointsx.push_back(posex[j].x_);
//      pointsy.push_back(posey[j].x_);
//  }
    int pnum =0,pnumx=0;
//  cout<< " 1111"<<endl;
    for(int i=0;i<posex.chains.size();i++)
    {
//      cout<< "chainx size: "<<posex.chains.size()<<endl;
//      cout<< "chainy size: "<<posey.chains.size()<<endl;
        Chain Chanx = posex.chains[i];
//      Chain Chany = posey.chains[i];
        for(int j=0;j<Chanx.residues.size();j++)
        {
//          cout<< "residuex size: "<< Chanx.residues.size()<<endl;
//          cout<< "residuey size: "<< Chany.residues.size()<<endl;
            Residue Resdx = Chanx.residues[j];
    //      Residue Resdy = Chany.residues[j];
            Rot rotx;
//          cout<<"Res: "<<Resdx.atoms.size()<<endl;
            pointsrenm.push_back(Resdx.resname);
//          pointsrenmy.push_back(Resdy.resname);
            for(int t=0;t<Resdx.atoms.size();t++)
            {
//              cout<< "atomx size: "<<Resdx.atoms.size()<<endl;
                poseCoord pcdx;
                Atom Atmx=Resdx.atoms[t];
                if(Atmx.atona != " CA " && Atmx.atona != " C  " && Atmx.atona != " N  " && Atmx.atona != " O  ")
                {
                    rotx.res_atm.push_back(Atmx.xyzVect);
                    vector<string> atmnmx = string_splitxx(Atmx.atona,' ');
                    string atmnmy = atmnmx[0];
                    rotx.resatmnm.push_back(atmnmy);
                }               
                if(Atmx.atona == " CA ") {rotx.ca = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);pnumx = pnumx +1;}
                if(Atmx.atona == " C  ") {rotx.cc = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);}//pnumx = pnumx +1;
                if(Atmx.atona == " N  ") {rotx.cn = Atmx.xyzVect;pcdx.x_=Atmx.xyzVect;pcdx.B_=Atmx.tempfac;pcdx.elt_= name2eltx(Atmx.atona);pointsBx.push_back(pcdx);}//pnumx = pnumx +1;
            /*    if(Atmx.atona == " C  ") {rotx.cc = Atmx.xyzVect;}
                if(Atmx.atona == " N  ") {rotx.cn = Atmx.xyzVect;}
                if(Atmx.atona == " O  ") {rotx.co = Atmx.xyzVect;}
                if(Atmx.atona == " CB ") {rotx.cb = Atmx.xyzVect;}
                if(Atmx.atona == " CG ") {rotx.cg = Atmx.xyzVect;}
                if(Atmx.atona == " CG1") {rotx.cg1 = Atmx.xyzVect;}
                if(Atmx.atona == " CG2") {rotx.cg2 = Atmx.xyzVect;}
                if(Atmx.atona == " OG ") {rotx.og  = Atmx.xyzVect;}
                if(Atmx.atona == " SG ") {rotx.sg = Atmx.xyzVect;}
                if(Atmx.atona == " CD ") {rotx.cd = Atmx.xyzVect;}
                if(Atmx.atona == " SD ") {rotx.sd = Atmx.xyzVect;}
                if(Atmx.atona == " CD1") {rotx.cd1 = Atmx.xyzVect;}
                if(Atmx.atona == " CD2") {rotx.cd2 = Atmx.xyzVect;}
                if(Atmx.atona == " ND2") {rotx.nd2 = Atmx.xyzVect;}
                if(Atmx.atona == " ND1") {rotx.nd1 = Atmx.xyzVect;}
                if(Atmx.atona == " OD2") {rotx.od2 = Atmx.xyzVect;}
                if(Atmx.atona == " OE2") {rotx.oe2 = Atmx.xyzVect;}
                if(Atmx.atona == " CE ") {rotx.ce = Atmx.xyzVect;}
                if(Atmx.atona == " NE2") {rotx.ne2 = Atmx.xyzVect;}
                if(Atmx.atona == " NE ") {rotx.ne = Atmx.xyzVect;}
                if(Atmx.atona == " NZ ") {rotx.nz = Atmx.xyzVect;}
                if(Atmx.atona == " CZ ") {rotx.cz = Atmx.xyzVect;}             */     
    /*          if(Atmy.atona == " CB ") roty.cb = Atmy.xyzVect;                    
                if(Atmx.atona == " CG ") rotx.cg = Atmy.xyzVect;
                if(Atmy.atona == " CG ") roty.cg = Atmy.xyzVect;
                if(Atmx.atona == " CD ") rotx.cd = Atmy.xyzVect;
                if(Atmy.atona == " CD ") roty.cd = Atmy.xyzVect;
                if(Atmx.atona == " CE ") rotx.ce = Atmy.xyzVect;
                if(Atmy.atona == " CE ") roty.ce = Atmy.xyzVect;
                if(Atmx.atona == " CZ ") rotx.cz = Atmy.xyzVect;
                if(Atmy.atona == " CZ ") roty.cz = Atmy.xyzVect;        */                  

            }
    /*      for(int t=0;t<Resdy.atoms.size();t++)
            {
//              cout<< "atomy size: "<<Resdy.atoms.size()<<endl;
                Atom Atmy=Resdy.atoms[t];
                if(Atmy.atona != " CA " && Atmy.atona != " C  " && Atmy.atona != " N  " && Atmy.atona != " O  ")
                {
                    roty.res_atm.push_back(Atmy.xyzVect);
                }               
                if(Atmy.atona == " CA ") {roty.ca = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
                if(Atmy.atona == " C  ") {roty.cc = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
                if(Atmy.atona == " N  ") {roty.cn = Atmy.xyzVect;pointsABy.push_back(Atmy.xyzVect);}
                if(Atmy.atona == " O  ") {roty.co = Atmy.xyzVect;}
                if(Atmy.atona == " CB ") {roty.cb = Atmy.xyzVect;}
                if(Atmy.atona == " CG ") {roty.cg = Atmy.xyzVect;}
                if(Atmy.atona == " CG1") {roty.cg1 = Atmy.xyzVect;}
                if(Atmy.atona == " CG2") {roty.cg2 = Atmy.xyzVect;}
                if(Atmy.atona == " OG ") {roty.og  = Atmy.xyzVect;}
                if(Atmy.atona == " SG ") {roty.sg = Atmy.xyzVect;}
                if(Atmy.atona == " CD ") {roty.cd = Atmy.xyzVect;}
                if(Atmy.atona == " SD ") {roty.sd = Atmy.xyzVect;}
                if(Atmy.atona == " CD1") {roty.cd1 = Atmy.xyzVect;}
                if(Atmy.atona == " CD2") {roty.cd2 = Atmy.xyzVect;}
                if(Atmy.atona == " ND2") {roty.nd2 = Atmy.xyzVect;}
                if(Atmy.atona == " ND1") {roty.nd1 = Atmy.xyzVect;}
                if(Atmy.atona == " OD2") {roty.od2 = Atmy.xyzVect;}
                if(Atmy.atona == " OE2") {roty.oe2 = Atmy.xyzVect;}
                if(Atmy.atona == " CE ") {roty.ce = Atmy.xyzVect;}
                if(Atmy.atona == " NE2") {roty.ne2 = Atmy.xyzVect;}
                if(Atmy.atona == " NE ") {roty.ne = Atmy.xyzVect;}
                if(Atmy.atona == " NZ ") {roty.nz = Atmy.xyzVect;}
                if(Atmy.atona == " CZ ") {roty.cz = Atmy.xyzVect;}              
            }           */
            pointsx.push_back(rotx);
//          pointsy.push_back(roty);
            pnum = pnum +1; 
//          cout<< pnum<<endl;      
        }
    }
    cout<<"pnum,pnumx: "<<pnum<<" "<<pnumx<<endl;
    pointsBx0 = pointsBx;

    vector<float> axyz;
    axyz=vector<float>(3,0.0);
    float ang; 
    float angle_rotate; 
    float old_dE1=0.0;
    float new_dE1=0.0;    
    float old_CC=0.0;
    float new_CC=0.0;
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
//    vector<poseCoord> fin_mat;
//    vector<vector<float> > fin_matx;
    vector<float> coord_change(3,0.0);

    vector<vector<float>> atm_idx(pnum,vector<float>(3,0.0));


    float effReso = std::max( 2.4+0.8*MRC_R , double(MRC_R ));
    float k=(M_PI/effReso)*(M_PI/effReso);
    float a=33.0;  // treat everything as ALA
    float C=a*pow(k/3.1415926,1.5);
    int tmp_i=0;
    ObjexxFCL::FArray3D< float > rhoC0;
    ObjexxFCL::FArray3D< float > rhoC01;
    ObjexxFCL::FArray3D< float > inv_rho_mask0;
    ObjexxFCL::FArray3D< float > inv_rho_mask1;
    rhoC0.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    rhoC01.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    inv_rho_mask0.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    inv_rho_mask1.dimension(theDensityMap.density.u1() , theDensityMap.density.u2() , theDensityMap.density.u3());
    for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
        rhoC0[t]=0.0;
        inv_rho_mask0[t]=1.0;
    }
    vector<float> del_ijx(3,0.0);
    vector<float> atm_jx(3,0.0); 
    string elt_i;    
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
            k = sig_j.k( theDensityMap.effectiveB );
            C = sig_j.C( k );
            if ( C < 1e-6 ) continue;   

            cartX1 = pointsBx[i].x_;           
            MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

            // the location of atom in grid ?
            atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
            atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
            atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
    //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
    //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
    //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
    //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
            for(int z=1;z<=theDensityMap.density.u3();z++)
            {
                atm_jx[2] = z;
                del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMap.grid[2];
                if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                del_ijx[0] = del_ijx[1] = 0.0;
                vector<float> frac_tmpz;
                MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                for(int y=1;y<=theDensityMap.density.u2();y++)
                {
                    atm_jx[1] = y;
                    del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMap.grid[1];
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
                        del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMap.grid[0];
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
    float max_den = 0.0;
    float max_rhc = 0.0;
    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
//        cout<<theDensityMap.density[x]<<" ";
        if(theDensityMap.density[x]>max_den) max_den = theDensityMap.density[x] ;
        if(rhoC0[x]>max_rhc) max_rhc =rhoC0[x] ;
    }
/*    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
//        cout<<theDensityMap.density[x]<<" ";
        if(theDensityMap.density[x]<(3.0/5.0)*max_den) theDensityMap.density[x]=0; 
    }    */
/*    for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
//        cout<<theDensityMap.density[x]<<" ";
        if(rhoC0[x]<(3.0/5.0)*max_rhc) rhoC0[x]=0; 
    }     */

   
    // ss information
/*    vector<string> ss_inf0(pnum);
    vector<vector<float> > CA_cor0(pnum,vector<float>(3,0.0));
    cout<<"RRRRRzz"<<endl;
    tmp_i=0;
    for(int i=0; i<pointsBx.size();i++)
    {
//        cout<<CA_cor0[tmp_i][0]<<'\t'<<CA_cor0[tmp_i][1]<<'\t'<<CA_cor0[tmp_i][2]<<endl;
        if(pointsBx[i].elt_ == "CA")
        {
            CA_cor0[tmp_i]=pointsBx[i].x_;
//            cout<<CA_cor0[tmp_i][0]<<'\t'<<CA_cor0[tmp_i][1]<<'\t'<<CA_cor0[tmp_i][2]<<endl;
            tmp_i=tmp_i+1;
//            cout<<pointsBx[i].x_[0]<<'\t'<<pointsBx[i].x_[1]<<'\t'<<pointsBx[i].x_[2]<<endl;
        }
    }
    make_secx(CA_cor0,pnum,ss_inf0);  */
    vector<boneinfo> bb(pnum);
    vector<poseCoord> proseq;
    proseq=pointsBx;
    int numbb=pnum;
    int tmp_tt=0;
    for(int j=0;j<pnum;j++)
    {
        bb[tmp_tt].indn=3*tmp_tt+0;
        bb[tmp_tt].indca=3*tmp_tt+1;
        bb[tmp_tt].indc=3*tmp_tt+2;   
        tmp_tt=tmp_tt+1;
    }     
    cout<<"xxExxx"<<endl;
    calcsse2(bb, numbb, proseq);
    cout<<"xxxxx"<<endl;

    for(int i=0;i<pnum;i++)
    {
        cout<<numtoss(bb[i].sst)<<" ";
    }    
    cout<<endl;
/*    for(int i=0;i<pnum;i++)
    {
        if(i==0) continue;
        int tyx=i-1;
        int tyy=i;
        int tyz=i+1;
        if(tyz>=pnum-1) break;
        if(bb[tyy].sst!=bb.sst[tyx] && bb[tyy].sst!=bb[tyz].sst)  bb[i].sst='C';
    } */
    for(int i=0;i<pnum;i++)
    {
        cout<<i<<":"<<numtoss(bb[i].sst)<<" ";
    }    
    cout<<endl;    

    // domain align
    // rotate around point
    int beg_p=dm_beg;
    int end_p=dm_end;
    int randop=-1; // rotate point
    int randopx=-1;
    if(beg_p>3) randop=beg_p;
    if(end_p<(pnum-4) && end_p>4) randopx=end_p;
    float bond0max=0.0,bond1max=0.0,bond2max=0.0;
    float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;   
    cout<<"CCCC"<<endl; 

    vector<float> mass_pdb(3,0.0);
    vector<float> mass_mrc(3,0.0);
    vector<float> coor_pdb(3,0.0);
    float allpdbmass=0.0;
    float allmrcmass=0.0;
    float olddis_p_m=0.0;
    float newdis_p_m=0.0;
    int acc_rate=0; // not acc_rate;
    for(int tttt=0;tttt<2;tttt++)
    {
        if(tttt==0)
        {
           if(dm_beg<3)
           {
                beg_p=0;
                end_p=dm_end;
           }else
           {
                beg_p=0;
                end_p=dm_beg;
           }

//            if(beg_p>3) randop=beg_p;
//            if(end_p<(pnum-4) && end_p>4) randopx=end_p;            
        }
        if(tttt==1)
        {
           if(dm_beg<3)
           {
                beg_p=dm_end;
                end_p=pnum;
           }else
           {
                beg_p=dm_beg;
                end_p=pnum;
           }
//            if(beg_p>3) randop=beg_p;
//            if(end_p<(pnum-4) && end_p>4) randopx=end_p;            
        }                
        float Tend[5]={0.01,0.003,0.008,0.001};
    vector<poseCoord> best_model;
    float best_E=10.0;        
    for(int ttt=0;ttt<4;ttt++)
    {
        float KT0 =0.01;
        float KTn =0.001;
        float KTx= pow(float(KTn/KT0),float(float(ttt)/float(3)));
    //    KT=KT0*float(KTx);
        KT=Tend[ttt];
        cout<<"KT: "<<KT<<endl;
    //    KT=0.005;

//        float angle0 =180.0;
        float angle0 =30.0;
//        float anglen =10;
        float anglen =5.0;
        float anglex= pow(float(anglen/angle0),float(float(ttt)/float(3)));
        float angle0n=angle0*float(anglex);
    //    float angle0n=30.0;

        float trdis0 =3.0;
        float trdisn=0.5;
        float trdisx= pow(float(trdisn/trdis0),float(float(ttt)/float(3)));
        float trdis0n=trdis0*float(trdisx);
    //    float trdis0n=2.0;      

    //    cout<<"trdisx,anglex,KT 0: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;  

        float dis_p=0.0;
        bool tp_clashx= false;
        old_vwd =0.0;
        old_clash = 0.0;
        old_dE =0.0;
        new_dE =0.0;
        axyz=vector<float>(3,0.0);  
        dE=0.0;  
        bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
        bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;            
    //    vector<vector<float> > fin_matx;
        float new_CC1=0.0;
        float old_CC1=0.0;
        float new_CC2=0.0;
        new_dE2=0.0;
        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0); 
        tmp_i=0;     
        for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
            rhoC01[t]=0.0;
            inv_rho_mask1[t]=1.0;
        }   

        vector<poseCoord > fin_maty;
        coor_pdb=vector<float>(3,0.0);
        for(int j=beg_p;j<end_p;j++)
        {
            fin_maty.push_back(pointsBx[3*j+0]);
            fin_maty.push_back(pointsBx[3*j+1]);
            fin_maty.push_back(pointsBx[3*j+2]);
            coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
            coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
            coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];
        }
        coor_pdb[0] = coor_pdb[0]/(end_p-beg_p);
        coor_pdb[1] = coor_pdb[1]/(end_p-beg_p);
        coor_pdb[2] = coor_pdb[2]/(end_p-beg_p);
        cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl; 
        atm_idx=vector<vector<float> > (fin_maty.size()-1,vector<float>(3,0.0));
        for(int i=0;i<fin_maty.size();i++)  // rotate fragment from beg_p to end_p
        {
             if(fin_maty[i].elt_ == "CA")
            {           
    //            cout<<"RRRR: "<<i;
                vector<float> cartX1;
                vector<float> fracX1;
                elt_i = fin_maty[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMap.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_maty[i].x_;           
                MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMap.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMap.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMap.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMap.grid[1];
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
                            del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMap.grid[0];
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
                            rhoC01(x,y,z) += atm;
                            inv_rho_mask1(x,y,z) *= (1 - inv_msk);


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
        } 

        // bandlimit mask at 'radius'
        ObjexxFCL::FArray3D< std::complex<double> > Finv_rho_mask1;
        fourier::fft3(inv_rho_mask1, Finv_rho_mask1);
        //int H,K,L;
        for ( int z=1; z<=(int)theDensityMap.grid[2]; ++z ) {
            int H = (z < (int)theDensityMap.grid[2]/2) ? z-1 : z-theDensityMap.grid[2] - 1;
            for ( int y=1; y<=(int)theDensityMap.grid[1]; ++y ) {
                int K = (y < (int)theDensityMap.grid[1]/2) ? y-1 : y-theDensityMap.grid[1]-1;
                for ( int x=1; x<=(int)theDensityMap.grid[0]; ++x ) {
                    int L = (x < (int)theDensityMap.grid[0]/2) ? x-1 : x-theDensityMap.grid[0]-1;
                    float S2c =  theDensityMap.S2(H,K,L);

                    // exp fade
                    float highres_limit=1000.0;
                    float scale = exp(-S2c*(highres_limit*highres_limit));
                    Finv_rho_mask1(x,y,z) *= scale;
                }
            }
        }
        fourier::ifft3(Finv_rho_mask1, inv_rho_mask1);      
          
        float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
        float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
        float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
        allpdbmass=0.0;
        allmrcmass=0.0;
        mass_pdb=vector<float>(3,0.0);
        mass_mrc=vector<float>(3,0.0);
        for(int x=int(bond0min);x<=bond0max;x++)
        {
           for(int y=int(bond1min);y<=bond1max;y++)
            {
                for(int z=int(bond2min);z<=bond2max;z++)
                {
                    clc_x2 = rhoC01(x,y,z);
    //                clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMap.density(x,y,z);
    //                obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask1(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                    mass_pdb[0]=  mass_pdb[0] + eps_x2*clc_x2*x;
                    mass_pdb[1]=  mass_pdb[1] + eps_x2*clc_x2*y;
                    mass_pdb[2]=  mass_pdb[2] + eps_x2*clc_x2*z;
                    mass_mrc[0]=  mass_mrc[0] + eps_x2*obs_x2*x;
                    mass_mrc[1]=  mass_mrc[1] + eps_x2*obs_x2*y;
                    mass_mrc[2]=  mass_mrc[2] + eps_x2*obs_x2*z;
                    allpdbmass = allpdbmass + eps_x2*clc_x2;
                    allmrcmass = allmrcmass + eps_x2*obs_x2;
                }
            }
        }
        mass_pdb[0] = mass_pdb[0]/allpdbmass;
        mass_pdb[1] = mass_pdb[1]/allpdbmass;
        mass_pdb[2] = mass_pdb[2]/allpdbmass;
        mass_mrc[0] = mass_mrc[0]/allmrcmass;
        mass_mrc[1] = mass_mrc[1]/allmrcmass;
        mass_mrc[2] = mass_mrc[2]/allmrcmass;
    //    olddis_p_m=Distance_point(mass_pdb,mass_mrc);   
        /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
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
                    }   */           
        varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
        varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 );
        if ( varC_i2 == 0 || varO_i2 == 0 ) {
            CC_i2 = 0;
        } else {
            CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
        }   
        new_CC2=1.0-CC_i2;   
    //    acc_rate=0;

        for(int iii=0;iii<500;iii++)
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
            float t0=trdis0n;
        //    float t0=1.0;

        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

            // Rotation matrix
            float anggg=angle0n;
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
            
            float asin=sin(angle_rotategg);
            float acos=cos(angle_rotategg);
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
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];

        //    beg_p=rand_point;   

            int tmp3=0;
            vector<poseCoord > fin_mat;
        //    fin_mat = tmp_mat;
        //    if(beg_p==0)
        //    {
                tmp3=0;
                vector<poseCoord>().swap(fin_mat);
                for(int j=beg_p;j<end_p;j++)
                {
                    fin_mat.push_back(pointsBx[3*j+0]);
                    fin_mat.push_back(pointsBx[3*j+1]);
                    fin_mat.push_back(pointsBx[3*j+2]);
                    tmp3=tmp3+1;
                }
                dis_p=0.0;
                tp_clashx= false;
                old_vwd =0.0;
                old_clash = 0.0;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                            //    if(abs(jj+3*beg_p-t)<1) continue;
                            //    if(abs(t)<3*end_p) continue;
                                if(abs(t)>(3*beg_p-1) && abs(t)<3*(end_p)) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    old_clash = old_clash + 1.0/ sqrt(dis_p) ;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }    
                float d1=0.0;
                float d2=0.0;
                if(randop>=0)
                {
                    d1= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[1].x_); // the distance between the beg_p CA atom and (beg_p-1) CA
                    d2= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[4].x_); // the distance between the beg_p CA atom and (beg_p-2) CA
                //    if(d2<d1) continue;
                } 
                if(randopx>=0)
                {
                    d1= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-2].x_); // the distance between the end_p CA atom and (end-1) CA
                    d2= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-5].x_); // the distance between the end_p CA atom and (end-2) CA
                //    if(d2<d1) continue;
                }                 
                old_CC1=new_CC2;   
                old_dE = old_CC1 ;// + old_clash + 0.1*d1/d2;// + olddis_p_m;
            //    old_dE = old_clash + olddis_p_m;                      
                vector<poseCoord > fin_matx=fin_mat;
                tmp3=0;
                for(int j=0;j<fin_mat.size();j++) // not last N
                {
                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                    tmp3=tmp3+1;
                }
                vector<poseCoord>().swap(fin_matx);
                dis_p=0.0;
                new_vwd =0.0;
                new_clash = 0.0;
            //    cout<<"size: "<<fin_mat.size()<<" "<<pointsBx.size()<<endl;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                            //    if(abs(t)<3*end_p) continue;
                                if(abs(t)>(3*beg_p-1) && abs(t)<3*(end_p)) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    new_clash = new_clash + 1.0/ sqrt(dis_p) ;
                 //                   cout<<"jj,t: "<<jj<<" "<<t<<" "<<dis_p<<endl;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }  

        if(randop>=0)
        {
            d1= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[1].x_); // the distance between the beg_p CA atom and (beg_p-1) CA
            d2= Distance_point(pointsBx[3*(randop-1)+1].x_,fin_mat[4].x_); // the distance between the beg_p CA atom and (beg_p-2) CA
        //    if(d2<d1) continue;
        } 
        if(randopx>=0)
        {
            d1= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-2].x_); // the distance between the end_p CA atom and (end-1) CA
            d2= Distance_point(pointsBx[3*(randopx)+1].x_,fin_mat[fin_mat.size()-5].x_); // the distance between the end_p CA atom and (end-2) CA
        //    if(d2<d1) continue;
        }     

        //    } 
        /*    else 
            {
                tmp3=0;
                vector<poseCoord>().swap(fin_mat);
                for(int j=beg_p;j<end_p;j++)
                {
                    fin_mat.push_back(pointsBx[3*j+0]);
                    fin_mat.push_back(pointsBx[3*j+1]);
                    fin_mat.push_back(pointsBx[3*j+2]);
                    tmp3=tmp3+1;
                } 
            //    cout<<"fin_mat coor0: "<<fin_mat[1].x_[0]<<" "<<fin_mat[1].x_[1]<<" "<<fin_mat[1].x_[2]<<endl;
                dis_p=0.0;
                old_vwd =0.0;
                old_clash = 0.0;
        //        cout<<"size: "<<fin_mat.size()<<" "<<pointsBx.size()<<endl;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                                if(abs(t)>=3*beg_p) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    old_clash = old_clash + 3.75*1.0/ sqrt(dis_p) ;
                            //        cout<<"jj,t,dis_p: "<<jj<<" "<<t<<" "<<dis_p<<endl;
                                //    cout<<"jjd: "<<fin_mat[jj].x_[0]<<" "<<fin_mat[jj].x_[1]<<" "<<fin_mat[jj].x_[2]<<endl;
                                //    cout<<"pd: "<<pointsBx[t].x_[0]<<" "<<pointsBx[t].x_[1]<<" "<<pointsBx[t].x_[2]<<endl;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }    
                old_CC1=new_CC2;   
                old_dE = old_CC1 + old_clash;//+ olddis_p_m; 
            //     old_dE =  old_clash+ olddis_p_m;                   
                vector<poseCoord > fin_matx=fin_mat;    
                tmp3=0;        
                for(int j=1;j<fin_mat.size();j++) // not first N
                {
                    fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                    fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                    fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                    tmp3=tmp3+1;
                }
                dis_p=0.0;
                new_vwd =0.0;
                new_clash = 0.0;
        //        cout<<"size: "<<fin_mat.size()<<" "<<pointsBx.size()<<endl;
                for(int jj=0;jj<fin_mat.size();jj++)
                {
                    if(fin_mat[jj].elt_=="CA")
                    {
                        for(int t=0;t<pointsBx.size();t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                                if(abs(t)>=3*beg_p) continue;
                                dis_p = Distance_point(fin_mat[jj].x_,pointsBx[t].x_);
                                if(dis_p < 3.70)
                                {
                                    new_clash = new_clash + 3.75*1.0/ sqrt(dis_p) ;
                 //                   cout<<"jj,t: "<<jj<<" "<<t<<" "<<dis_p<<endl;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                            }
                        }
                    }
                }                                         
            }       */

        //    cout<<"fin_mat coor1: "<<fin_mat[1].x_[0]<<" "<<fin_mat[1].x_[1]<<" "<<fin_mat[1].x_[2]<<endl;  
        //    if(new_clash>10) continue;
            bond0max=-10000.0;bond1max=-10000.0;bond2max=-10000.0;
            bond0min=10000.0;bond1min=10000.0;bond2min=10000.0;            
        //    vector<vector<float> > fin_matx;
            del_ijx=vector<float>(3,0.0);
            atm_jx=vector<float>(3,0.0); 
            tmp_i=0;     
            for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
                rhoC01[t]=0.0;
                inv_rho_mask1[t]=1.0;
            }   
            atm_idx=vector<vector<float> > (fin_mat.size()-1,vector<float>(3,0.0));
            for(int i=0;i<fin_mat.size();i++)  // rotate fragment from beg_p to end_p
            {
                 if(fin_mat[i].elt_ == "CA")
                {           
        //            cout<<"RRRR: "<<i;
                    vector<float> cartX1;
                    vector<float> fracX1;
                    elt_i = fin_mat[i].elt_;
                    elt_i = elt_i[0];
                    OneGaussianScattering sig_j = get_A( elt_i );
                    k = sig_j.k( theDensityMap.effectiveB );
                    C = sig_j.C( k );
                    if ( C < 1e-6 ) continue;   

                    cartX1 = fin_mat[i].x_;           
                    MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                    // the location of atom in grid ?
                    atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
            //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
            //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
            //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
            //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                    for(int z=1;z<=theDensityMap.density.u3();z++)
                    {
                        atm_jx[2] = z;
                        del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMap.grid[2];
                        if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                        if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                        del_ijx[0] = del_ijx[1] = 0.0;
                        vector<float> frac_tmpz;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                        if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                        for(int y=1;y<=theDensityMap.density.u2();y++)
                        {
                            atm_jx[1] = y;
                            del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMap.grid[1];
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
                                del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMap.grid[0];
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
                                rhoC01(x,y,z) += atm;
                                inv_rho_mask1(x,y,z) *= (1 - inv_msk);


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
            } 

            ObjexxFCL::FArray3D< std::complex<double> > Finv_rho_mask2;
            fourier::fft3(inv_rho_mask1, Finv_rho_mask2);
            //int H,K,L;
            for ( int z=1; z<=(int)theDensityMap.grid[2]; ++z ) {
                int H = (z < (int)theDensityMap.grid[2]/2) ? z-1 : z-theDensityMap.grid[2] - 1;
                for ( int y=1; y<=(int)theDensityMap.grid[1]; ++y ) {
                    int K = (y < (int)theDensityMap.grid[1]/2) ? y-1 : y-theDensityMap.grid[1]-1;
                    for ( int x=1; x<=(int)theDensityMap.grid[0]; ++x ) {
                        int L = (x < (int)theDensityMap.grid[0]/2) ? x-1 : x-theDensityMap.grid[0]-1;
                        float S2c =  theDensityMap.S2(H,K,L);

                        // exp fade
                        float highres_limit=1000.0;
                        float scale = exp(-S2c*(highres_limit*highres_limit));
                        Finv_rho_mask2(x,y,z) *= scale;
                    }
                }
            }
            fourier::ifft3(Finv_rho_mask2, inv_rho_mask1);

            sumC_i2=0.0; sumO_i2=0.0; sumCO_i2=0.0; vol_i2=0.0; CC_i2=0.0;
            sumO2_i2=0.0; sumC2_i2=0.0; varC_i2=0.0; varO_i2=0.0;
            clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
            allpdbmass=0.0;
            allmrcmass=0.0;
            mass_pdb=vector<float>(3,0.0);
            mass_mrc=vector<float>(3,0.0);            
            for(int x=int(bond0min);x<=bond0max;x++)
            {
               for(int y=int(bond1min);y<=bond1max;y++)
                {
                    for(int z=int(bond2min);z<=bond2max;z++)
                    {
                        clc_x2 = rhoC01(x,y,z);
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density(x,y,z);
        //                obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask1(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                        sumCO_i2 += eps_x2*clc_x2*obs_x2;
                        sumO_i2  += eps_x2*obs_x2;
                        sumO2_i2 += eps_x2*obs_x2*obs_x2;
                        sumC_i2  += eps_x2*clc_x2;
                        sumC2_i2 += eps_x2*clc_x2*clc_x2;
                        vol_i2   += eps_x2;
                        mass_pdb[0]=  mass_pdb[0] + eps_x2*clc_x2*x;
                        mass_pdb[1]=  mass_pdb[1] + eps_x2*clc_x2*y;
                        mass_pdb[2]=  mass_pdb[2] + eps_x2*clc_x2*z;
                        mass_mrc[0]=  mass_mrc[0] + eps_x2*obs_x2*x;
                        mass_mrc[1]=  mass_mrc[1] + eps_x2*obs_x2*y;
                        mass_mrc[2]=  mass_mrc[2] + eps_x2*obs_x2*z;
                        allpdbmass = allpdbmass + eps_x2*clc_x2;
                        allmrcmass = allmrcmass + eps_x2*obs_x2;                        
                    }
                }
            }                       
            varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2 / vol_i2 );
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }   
            mass_pdb[0] = mass_pdb[0]/allpdbmass;
            mass_pdb[1] = mass_pdb[1]/allpdbmass;
            mass_pdb[2] = mass_pdb[2]/allpdbmass;
            mass_mrc[0] = mass_mrc[0]/allmrcmass;
            mass_mrc[1] = mass_mrc[1]/allmrcmass;
            mass_mrc[2] = mass_mrc[2]/allmrcmass;
//            newdis_p_m=Distance_point(mass_pdb,mass_mrc);             
/*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
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
            }   */           
        // new clash
        //    cout<<"new_clash: "<<new_clash<<" ";           
            new_CC1=1.0-CC_i2; 
            new_dE = new_CC1 ;// + new_clash+0.1*d1/d2;// + newdis_p_m;
        //    new_dE = new_clash + newdis_p_m;
       //    cout<<"old_CC1,new_CC1:R "<<old_CC1<<" "<<new_CC1<<endl;
            cout<<"old_CC1,new_CC1:RRRRRRRRRRRRR "<<old_CC1<<" "<<new_CC1<<endl;
            dE = new_dE-old_dE;
            acc_rate=acc_rate+1;
            if(new_dE<old_dE)
            {
                tmp3=0;
                coor_pdb=vector<float>(3,0.0);
                for(int j=beg_p;j<end_p;j++)
                {
                    pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                    pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                    pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                    coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                    coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                    coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                    
                    tmp3=tmp3+1;
                }
                coor_pdb[0] = coor_pdb[0]/(end_p-beg_p);
                coor_pdb[1] = coor_pdb[1]/(end_p-beg_p);
                coor_pdb[2] = coor_pdb[2]/(end_p-beg_p); 
                if(new_dE<best_E)
                {
                    best_E = new_dE;
                    vector<poseCoord>().swap(best_model);
                    best_model = fin_mat;
                }
            //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                  
                new_CC2 = new_CC1;
                acc_rate=0;
                olddis_p_m = newdis_p_m;
                 
                cout<<"old_CC1,new_CC1: "<<old_dE<<" "<<new_dE<<endl;         
            }else
            {
            //    float tmpx=rand()/double(RAND_MAX);
                float tmpx=randf0and1();
                float mc_v = exp(-dE/(KT));
                if(tmpx< mc_v) // CC must be >0
                {
                    tmp3=0;
                    coor_pdb=vector<float>(3,0.0);
                    for(int j=beg_p;j<end_p;j++)
                    {
                        pointsBx[3*j+0]=fin_mat[3*tmp3+0];
                        pointsBx[3*j+1]=fin_mat[3*tmp3+1];
                        pointsBx[3*j+2]=fin_mat[3*tmp3+2];

                        coor_pdb[0]= coor_pdb[0] + fin_mat[3*tmp3+1].x_[0];
                        coor_pdb[1]= coor_pdb[1] + fin_mat[3*tmp3+1].x_[1];
                        coor_pdb[2]= coor_pdb[2] + fin_mat[3*tmp3+1].x_[2];                         
                        tmp3=tmp3+1;
                    } 
                    coor_pdb[0] = coor_pdb[0]/(end_p-beg_p);
                    coor_pdb[1] = coor_pdb[1]/(end_p-beg_p);
                    coor_pdb[2] = coor_pdb[2]/(end_p-beg_p); 
                //    cout<<"coor: "<<coor_pdb[0]<<" "<<coor_pdb[1]<<" "<<coor_pdb[2]<<endl;                     
                    new_CC2 = new_CC1;
                    acc_rate=0;
                    olddis_p_m = newdis_p_m;
                    cout<<"old_CC1,new_CC1:XX "<<old_dE<<" "<<new_dE<<endl;                    
                }                
            }           
        }      
    }
        int tmp3=0;
        for(int j=beg_p;j<end_p;j++)
        {
            pointsBx[3*j+0]=best_model[3*tmp3+0];
            pointsBx[3*j+1]=best_model[3*tmp3+1];
            pointsBx[3*j+2]=best_model[3*tmp3+2];              
            tmp3=tmp3+1;
        }    
    }

    cout<<"CCCCCC"<<endl;

    for(int jjx=0;jjx<0;jjx++)
    {
        float new_CC1=0.0;
        float old_CC1=0.0;
        float new_CC2=0.0;
        KT =0.001;
        ang=1.0;
        int rand_a=0,rand_b=0;
    //    int randtx = rand()%pnum;
        int randtx = randIntCustom(0,pnum-1);
        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
        {
            int jmx=randtx,jnx=randtx;
            int nn=0;
            for(int j=0;j<pnum;j++)
            {
                jnx=randtx+j;
                if(jnx>=(pnum-1)) break;
                if(jnx==(end_p-1)||jnx==(beg_p-1))
                {
                    break;
                }
                if(nn>2 && numtoss(bb[jnx].sst)!='C') break;
                if(numtoss(bb[jnx].sst)=='C') 
                {
                    nn=nn+1;
                    if(nn>=4) break;
                }
            }
            nn=0;
            for(int j=0;j<pnum;j++)
            {
                jmx=randtx-j;
                if(jmx<=0) break;
                if(jmx==(end_p+1)||jmx==(beg_p))
                {
                    break;
                }
                if(nn>2 && numtoss(bb[jmx].sst)!='C') break;
                if(numtoss(bb[jmx].sst)=='C') 
                {
                    nn=nn+1;
                    if(nn>=4) break;
                }
            }
            rand_a=jmx;
            rand_b=jnx;            
        } else
        {
        //     float randt = rand()%4+3;
            float randt = randIntCustom(0,3)+3;
            int jmx=randtx,jnx=randtx;
            for(int j=0;j<randt;j++)
            {
                jnx=randtx+j;
                if(jnx>=(pnum-1)) break;
        //        if(randtx<randop)
        //        {
        //            if(jnx>=(randop-1)) break;
        //        }                
                if(numtoss(bb[jnx].sst)!='C') break;
            }
            if((jnx-randtx+1)<randt)
            {
                for(int j=0;j<(randt-(jnx-randtx));j++)
                {
                    jmx=randtx-j;
                    if(jmx<=0) break;
        //            if(randtx>=randop)
        //            {
        //                if(jmx<=randop) break;
        //            }                    
                    if(numtoss(bb[jmx].sst)!='C') break;
                }
            }
            rand_a=jmx;
            rand_b=jnx;                
        }

        if(abs(rand_b-rand_a+1)<3) continue;
        cout<<"rand a b: "<<rand_a<<" "<<rand_b<<endl;
        vector<poseCoord > fin_maty;
        vector<poseCoord>().swap(fin_maty);
        coor_pdb=vector<float>(3,0.0);
        for(int j=rand_a;j<=rand_b;j++)
        {
            fin_maty.push_back(pointsBx[3*j+0]);
            fin_maty.push_back(pointsBx[3*j+1]);
            fin_maty.push_back(pointsBx[3*j+2]);
            coor_pdb[0]= coor_pdb[0] + pointsBx[3*j+1].x_[0];
            coor_pdb[1]= coor_pdb[1] + pointsBx[3*j+1].x_[1];
            coor_pdb[2]= coor_pdb[2] + pointsBx[3*j+1].x_[2];             
        }
        coor_pdb[0] = coor_pdb[0]/abs(rand_b-rand_a+1);
        coor_pdb[1] = coor_pdb[1]/abs(rand_b-rand_a+1);
        coor_pdb[2] = coor_pdb[2]/abs(rand_b-rand_a+1);    
        del_ijx=vector<float>(3,0.0);
        atm_jx=vector<float>(3,0.0); 
        tmp_i=0;     
        for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
            rhoC01[t]=0.0;
            inv_rho_mask1[t]=1.0;
        }
        atm_idx=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
        bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
        bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;        
        for(int i=0;i<fin_maty.size();i++)  // rotate fragment from beg_p to end_p
        {
             if(fin_maty[i].elt_ == "CA")
            {           
    //            cout<<"RRRR: "<<i;
                vector<float> cartX1;
                vector<float> fracX1;
                elt_i = fin_maty[i].elt_;
                elt_i = elt_i[0];
                OneGaussianScattering sig_j = get_A( elt_i );
                k = sig_j.k( theDensityMap.effectiveB );
                C = sig_j.C( k );
                if ( C < 1e-6 ) continue;   

                cartX1 = fin_maty[i].x_;           
                MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);

                // the location of atom in grid ?
                atm_idx[tmp_i][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                atm_idx[tmp_i][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                atm_idx[tmp_i][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);   
        //        atm_idx[tmp_i][0] = fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1 ;
        //        atm_idx[tmp_i][1] = fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1 ;
        //        atm_idx[tmp_i][2] = fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1 ;
        //        cout<<"RRRR: "<<atm_idx[tmp_i][0]<<" "<<atm_idx[tmp_i][1]<<" "<<atm_idx[tmp_i][2]<<endl; 
                for(int z=1;z<=theDensityMap.density.u3();z++)
                {
                    atm_jx[2] = z;
                    del_ijx[2] =(atm_idx[tmp_i][2]-atm_jx[2])/theDensityMap.grid[2];
                    if ( del_ijx[2] > 0.5 ) del_ijx[2]-=1.0;
                    if ( del_ijx[2] < -0.5 ) del_ijx[2]+=1.0;                    
                    del_ijx[0] = del_ijx[1] = 0.0;
                    vector<float> frac_tmpz;
                    MatrixTimesTransVector(theDensityMap.f2c,del_ijx,frac_tmpz);
                    if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.CA_MASK) ) continue;                    
                    for(int y=1;y<=theDensityMap.density.u2();y++)
                    {
                        atm_jx[1] = y;
                        del_ijx[1] = (atm_idx[tmp_i][1] - atm_jx[1])/theDensityMap.grid[1];
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
                            del_ijx[0] = (atm_idx[tmp_i][0] - atm_jx[0])/theDensityMap.grid[0];
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
                            rhoC01(x,y,z) += atm;
                            inv_rho_mask1(x,y,z) *= (1 - inv_msk);


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
        } 
        float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
        float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
        float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
//        allpdbmass=0.0;
//        allmrcmass=0.0;
//        mass_pdb=vector<float>(3,0.0);
//        mass_mrc=vector<float>(3,0.0);
        for(int x=int(bond0min);x<=bond0max;x++)
        {
           for(int y=int(bond1min);y<=bond1max;y++)
            {
                for(int z=int(bond2min);z<=bond2max;z++)
                {
                    clc_x2 = rhoC01(x,y,z);
    //                clc_x2 = tmp_rhc;
                    obs_x2 = theDensityMap.density(x,y,z);
    //                obs_x2 = tmp_den;
                    eps_x2 = 1-inv_rho_mask1(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
    //                mass_pdb[0]=  mass_pdb[0] + eps_x2*clc_x2*x;
    //                mass_pdb[1]=  mass_pdb[1] + eps_x2*clc_x2*y;
    //                mass_pdb[2]=  mass_pdb[2] + eps_x2*clc_x2*z;
    //                mass_mrc[0]=  mass_mrc[0] + eps_x2*obs_x2*x;
    //                mass_mrc[1]=  mass_mrc[1] + eps_x2*obs_x2*y;
    //                mass_mrc[2]=  mass_mrc[2] + eps_x2*obs_x2*z;
    //                allpdbmass = allpdbmass + eps_x2*clc_x2;
    //                allmrcmass = allmrcmass + eps_x2*obs_x2;
                }
            }
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

        float dis_p=0.0;
        old_vwd =0.0;
        old_clash = 0.0;
        for(int jj=0;jj<fin_maty.size();jj++)
        {
            if(fin_maty[jj].elt_=="CA")
            {
                for(int t=0;t<pointsBx.size();t++)
                {
                    if(pointsBx[t].elt_=="CA")
                    {
                        if(abs(jj+3*rand_a-t)<2) continue;
                    //    if(abs(t)<3*end_p) continue;
                        dis_p = Distance_point(fin_maty[jj].x_,pointsBx[t].x_);
                        if(dis_p < 3.70)
                        {
                            old_clash = old_clash + 1.0/ sqrt(dis_p) ;
                        //    cout<<"jj t: "<<jj<<" "<<t<<" "<<old_clash<<endl;
                        //    cout<<"coor0: "<<fin_maty[jj].x_[0]<<" "<<fin_maty[jj].x_[1]<<" "<<fin_maty[jj].x_[2]<<endl;
                        //    cout<<"coor1: "<<pointsBx[t].x_[0]<<" "<<pointsBx[t].x_[1]<<" "<<pointsBx[t].x_[2]<<endl;
        //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                        }
                    }
                }
            }
        }     
    //    old_dE = old_CC1 + old_clash;// + olddis_p_m;        
        old_dE = 0.5*old_CC1 ;//+ 0.5*old_clash/abs(rand_b-rand_a+1);

     /*   int rotx_point=0; // rotate point
        if(rand_a==0||rand_a==randop||rand_b==(randop-1)||rand_b==(pnum-1))
        {
            if(rand_a==0||rand_a==randop)
            {
                rotx_point = rand_b;
            }else
            {
                rotx_point = rand_a;
            }

            float asin_theta=2.0*RandomDoubleX(0.0,1.0)-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*RandomDoubleX(0.0,1.0);
            float awx=acos_theta*cos(apha);
            float awy=acos_theta*sin(apha);
            float awz=asin_theta;
            // Translation Vector
            float t0=0.1;
            float t1=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;
            float t2=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;
            float t3=(RandomDoubleX(0.0,1.0)*2.0-1.0)*t0+0.0;
        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

            // Rotation matrix
            float anggg=1.0;
            float angle_rotategg=(2.0*RandomDoubleX(0.0,1.0)-1.0)*anggg; // rotate angle
            float asin=sin(angle_rotategg);
            float acos=cos(angle_rotategg);
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

            int rand_point = rotx_point;
            axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
            axyz[1]=pointsBx[3*rand_point+1].x_[1];
            axyz[2]=pointsBx[3*rand_point+1].x_[2];     

            int tmp3=0;
            vector<poseCoord > fin_mat;
            vector<poseCoord>().swap(fin_mat);
            fin_mat = fin_maty;
            if(rand_a==0||rand_a==randop)
            {
                for(int j=0;j<fin_mat.size()-1;j++) // not last N
                {
                    fin_mat[tmp3].x_[0]=t1+axyz[0]+(fin_maty[j].x_[0]-axyz[0])*u[0][0]+(fin_maty[j].x_[1]-axyz[1])*u[0][1]+(fin_maty[j].x_[2]-axyz[2])*u[0][2];
                    fin_mat[tmp3].x_[1]=t2+axyz[1]+(fin_maty[j].x_[0]-axyz[0])*u[1][0]+(fin_maty[j].x_[1]-axyz[1])*u[1][1]+(fin_maty[j].x_[2]-axyz[2])*u[1][2];
                    fin_mat[tmp3].x_[2]=t3+axyz[2]+(fin_maty[j].x_[0]-axyz[0])*u[2][0]+(fin_maty[j].x_[1]-axyz[1])*u[2][1]+(fin_maty[j].x_[2]-axyz[2])*u[2][2];
                    tmp3=tmp3+1;
                }                
            } else
            {
                for(int j=1;j<fin_mat.size();j++) // not first N
                {
                    fin_mat[tmp3].x_[0]=t1+axyz[0]+(fin_maty[j].x_[0]-axyz[0])*u[0][0]+(fin_maty[j].x_[1]-axyz[1])*u[0][1]+(fin_maty[j].x_[2]-axyz[2])*u[0][2];
                    fin_mat[tmp3].x_[1]=t2+axyz[1]+(fin_maty[j].x_[0]-axyz[0])*u[1][0]+(fin_maty[j].x_[1]-axyz[1])*u[1][1]+(fin_maty[j].x_[2]-axyz[2])*u[1][2];
                    fin_mat[tmp3].x_[2]=t3+axyz[2]+(fin_maty[j].x_[0]-axyz[0])*u[2][0]+(fin_maty[j].x_[1]-axyz[1])*u[2][1]+(fin_maty[j].x_[2]-axyz[2])*u[2][2];
                    tmp3=tmp3+1;
                }                
            }           
            vector<poseCoord>().swap(fin_maty);
            fin_maty=fin_mat;                

        }
        else  */
        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E') // rotate around mass center
        {
    //        cout<<"trdisx,anglex,KT: "<<trdis0n<<" "<<angle0n<<" "<<KT<<endl;
            float asin_theta=2.0*randf0and1()-1.0;
            float acos_theta=sqrt(1.0-asin_theta*asin_theta);
            float apha=2.0*PI*randf0and1();
            float awx=acos_theta*cos(apha);
            float awy=acos_theta*sin(apha);
            float awz=asin_theta;
            // Translation Vector
            float t0=0.1;
            float t1=(randf0and1()*2.0-1.0)*t0+0.0;
            float t2=(randf0and1()*2.0-1.0)*t0+0.0;
            float t3=(randf0and1()*2.0-1.0)*t0+0.0;
        //    cout<<"t0: "<<t1<<" "<<t2<<" "<<t3<<endl;

            // Rotation matrix
            float anggg=1.0;
            float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
            float asin=sin(angle_rotategg);
            float acos=cos(angle_rotategg);
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
            axyz[0]= coor_pdb[0];  
            axyz[1]= coor_pdb[1];
            axyz[2]= coor_pdb[2];
            cout<<"cent: "<< randtx<<" "<< axyz[0]<<" "<<axyz[1]<<" "<< axyz[2]<<endl; 

        //    beg_p=rand_point;   

            int tmp3=0;
            vector<poseCoord > fin_mat;
        //    fin_mat = tmp_mat;
            tmp3=0;
            vector<poseCoord>().swap(fin_mat);
            for(int j=rand_a+2;j<=rand_b-2;j++)
            {
                fin_mat.push_back(pointsBx[3*j+0]);
                fin_mat.push_back(pointsBx[3*j+1]);
                fin_mat.push_back(pointsBx[3*j+2]);               
                tmp3=tmp3+1;
            }        
        //    old_dE = old_clash + olddis_p_m;                      
            vector<poseCoord > fin_matx=fin_mat;
            tmp3=0;
            for(int j=0;j<fin_mat.size();j++) // not last N
            {
                fin_mat[j].x_[0]=t1+axyz[0]+(fin_matx[j].x_[0]-axyz[0])*u[0][0]+(fin_matx[j].x_[1]-axyz[1])*u[0][1]+(fin_matx[j].x_[2]-axyz[2])*u[0][2];
                fin_mat[j].x_[1]=t2+axyz[1]+(fin_matx[j].x_[0]-axyz[0])*u[1][0]+(fin_matx[j].x_[1]-axyz[1])*u[1][1]+(fin_matx[j].x_[2]-axyz[2])*u[1][2];
                fin_mat[j].x_[2]=t3+axyz[2]+(fin_matx[j].x_[0]-axyz[0])*u[2][0]+(fin_matx[j].x_[1]-axyz[1])*u[2][1]+(fin_matx[j].x_[2]-axyz[2])*u[2][2];
                tmp3=tmp3+1;
            }


            int numseq=abs(rand_b-rand_a+1);
            if((rand_b-rand_a-3)<2) continue;
            vector<point3f> decstr(5);

                tmp_tt=0;
                int tmp_ttx=0;
                for(int j=rand_a;j<=rand_a+4;j++)
                {
                    if(tmp_tt<2)
                    {
                        decstr[tmp_tt].ss2=numtoss(bb[j].sst);
                        decstr[tmp_tt].x=fin_maty[3*tmp_tt+1].x_[0];
                        decstr[tmp_tt].y=fin_maty[3*tmp_tt+1].x_[1];
                        decstr[tmp_tt].z=fin_maty[3*tmp_tt+1].x_[2];
                        decstr[tmp_tt].ptn.x=fin_maty[3*tmp_tt+0].x_[0];
                        decstr[tmp_tt].ptn.y=fin_maty[3*tmp_tt+0].x_[1];
                        decstr[tmp_tt].ptn.z=fin_maty[3*tmp_tt+0].x_[2];
                        decstr[tmp_tt].ptc.x=fin_maty[3*tmp_tt+2].x_[0];
                        decstr[tmp_tt].ptc.y=fin_maty[3*tmp_tt+2].x_[1];
                        decstr[tmp_tt].ptc.z=fin_maty[3*tmp_tt+2].x_[2];
                    }else 
                    {
                        decstr[tmp_tt].ss2=numtoss(bb[j].sst);
                        decstr[tmp_tt].x=fin_mat[3*tmp_ttx+1].x_[0];
                        decstr[tmp_tt].y=fin_mat[3*tmp_ttx+1].x_[1];
                        decstr[tmp_tt].z=fin_mat[3*tmp_ttx+1].x_[2];
                        decstr[tmp_tt].ptn.x=fin_mat[3*tmp_ttx+0].x_[0];
                        decstr[tmp_tt].ptn.y=fin_mat[3*tmp_ttx+0].x_[1];
                        decstr[tmp_tt].ptn.z=fin_mat[3*tmp_ttx+0].x_[2];
                        decstr[tmp_tt].ptc.x=fin_mat[3*tmp_ttx+2].x_[0];
                        decstr[tmp_tt].ptc.y=fin_mat[3*tmp_ttx+2].x_[1];
                        decstr[tmp_tt].ptc.z=fin_mat[3*tmp_ttx+2].x_[2]; 
                        tmp_ttx=tmp_ttx+1;                       
                    }
                    tmp_tt=tmp_tt+1;
                }

            mcfragsweepLMP2(decstr,numseq,1,4);
                tmp_tt=0;
                tmp_ttx=0;
                for(int j=rand_a;j<=rand_a+4;j++)
                {
                    if(tmp_tt<2)
                    {
                        fin_maty[3*tmp_tt+1].x_[0]=decstr[tmp_tt].x;
                        fin_maty[3*tmp_tt+1].x_[1]=decstr[tmp_tt].y;
                        fin_maty[3*tmp_tt+1].x_[2]=decstr[tmp_tt].z;
                        fin_maty[3*tmp_tt+0].x_[0]=decstr[tmp_tt].ptn.x;
                        fin_maty[3*tmp_tt+0].x_[1]=decstr[tmp_tt].ptn.y;
                        fin_maty[3*tmp_tt+0].x_[2]=decstr[tmp_tt].ptn.z;
                        fin_maty[3*tmp_tt+2].x_[0]=decstr[tmp_tt].ptc.x;
                        fin_maty[3*tmp_tt+2].x_[1]=decstr[tmp_tt].ptc.y;
                        fin_maty[3*tmp_tt+2].x_[2]=decstr[tmp_tt].ptc.z;
                    }else 
                    {
                        fin_mat[3*tmp_ttx+1].x_[0]=decstr[tmp_tt].x;
                        fin_mat[3*tmp_ttx+1].x_[1]=decstr[tmp_tt].y;
                        fin_mat[3*tmp_ttx+1].x_[2]=decstr[tmp_tt].z;
                        fin_mat[3*tmp_ttx+0].x_[0]=decstr[tmp_tt].ptn.x;
                        fin_mat[3*tmp_ttx+0].x_[1]=decstr[tmp_tt].ptn.y;
                        fin_mat[3*tmp_ttx+0].x_[2]=decstr[tmp_tt].ptn.z;
                        fin_mat[3*tmp_ttx+2].x_[0]=decstr[tmp_tt].ptc.x;
                        fin_mat[3*tmp_ttx+2].x_[1]=decstr[tmp_tt].ptc.y;
                        fin_mat[3*tmp_ttx+2].x_[2]=decstr[tmp_tt].ptc.z; 
                        tmp_ttx=tmp_ttx+1;                       
                    }
                    tmp_tt=tmp_tt+1;
                }
                
            vector<point3f>().swap(decstr);
            vector<point3f> decstrx(5);
            tmp_tt=0;
            tmp_ttx=2;
            int tmp_tttx=3;
            int tmr_numx=abs(rand_b-rand_a-3);
            if((rand_b-rand_a-3)<0) continue;
            int tmr_num = abs(rand_b-rand_a+1);
                for(int j=rand_b-4;j<=rand_b;j++)
                {
                    if(tmp_tt>2)
                    {
                        decstrx[tmp_tt].ss2=numtoss(bb[j].sst);
                        decstrx[tmp_tt].x=fin_maty[3*(tmr_num-tmp_ttx)+1].x_[0];
                        decstrx[tmp_tt].y=fin_maty[3*(tmr_num-tmp_ttx)+1].x_[1];
                        decstrx[tmp_tt].z=fin_maty[3*(tmr_num-tmp_ttx)+1].x_[2];
                        decstrx[tmp_tt].ptn.x=fin_maty[3*(tmr_num-tmp_ttx)+0].x_[0];
                        decstrx[tmp_tt].ptn.y=fin_maty[3*(tmr_num-tmp_ttx)+0].x_[1];
                        decstrx[tmp_tt].ptn.z=fin_maty[3*(tmr_num-tmp_ttx)+0].x_[2];
                        decstrx[tmp_tt].ptc.x=fin_maty[3*(tmr_num-tmp_ttx)+2].x_[0];
                        decstrx[tmp_tt].ptc.y=fin_maty[3*(tmr_num-tmp_ttx)+2].x_[1];
                        decstrx[tmp_tt].ptc.z=fin_maty[3*(tmr_num-tmp_ttx)+2].x_[2];
                        tmp_ttx = tmp_ttx -1;
                    }else 
                    {
                        decstrx[tmp_tt].ss2=numtoss(bb[j].sst);
                        decstrx[tmp_tt].x=fin_mat[3*(tmr_numx-tmp_tttx)+1].x_[0];
                        decstrx[tmp_tt].y=fin_mat[3*(tmr_numx-tmp_tttx)+1].x_[1];
                        decstrx[tmp_tt].z=fin_mat[3*(tmr_numx-tmp_tttx)+1].x_[2];
                        decstrx[tmp_tt].ptn.x=fin_mat[3*(tmr_numx-tmp_tttx)+0].x_[0];
                        decstrx[tmp_tt].ptn.y=fin_mat[3*(tmr_numx-tmp_tttx)+0].x_[1];
                        decstrx[tmp_tt].ptn.z=fin_mat[3*(tmr_numx-tmp_tttx)+0].x_[2];
                        decstrx[tmp_tt].ptc.x=fin_mat[3*(tmr_numx-tmp_tttx)+2].x_[0];
                        decstrx[tmp_tt].ptc.y=fin_mat[3*(tmr_numx-tmp_tttx)+2].x_[1];
                        decstrx[tmp_tt].ptc.z=fin_mat[3*(tmr_numx-tmp_tttx)+2].x_[2]; 
                        tmp_tttx=tmp_tttx-1;                       
                    }
                    tmp_tt=tmp_tt+1;
                }            
            mcfragsweepLMP2(decstrx,numseq,1,4);

            tmp_tt=0;
            tmp_ttx=2;
            tmp_tttx=3;
            tmr_numx=abs(rand_b-rand_a-3);
            tmr_num = abs(rand_b-rand_a+1);
                for(int j=rand_b-4;j<=rand_b;j++)
                {
                    if(tmp_tt>2)
                    {
                        fin_maty[3*(tmr_num-tmp_ttx)+1].x_[0]=decstrx[tmp_tt].x;
                        fin_maty[3*(tmr_num-tmp_ttx)+1].x_[1]=decstrx[tmp_tt].y;
                        fin_maty[3*(tmr_num-tmp_ttx)+1].x_[2]=decstrx[tmp_tt].z;
                        fin_maty[3*(tmr_num-tmp_ttx)+0].x_[0]=decstrx[tmp_tt].ptn.x;
                        fin_maty[3*(tmr_num-tmp_ttx)+0].x_[1]=decstrx[tmp_tt].ptn.y;
                        fin_maty[3*(tmr_num-tmp_ttx)+0].x_[2]=decstrx[tmp_tt].ptn.z;
                        fin_maty[3*(tmr_num-tmp_ttx)+2].x_[0]=decstrx[tmp_tt].ptc.x;
                        fin_maty[3*(tmr_num-tmp_ttx)+2].x_[1]=decstrx[tmp_tt].ptc.y;
                        fin_maty[3*(tmr_num-tmp_ttx)+2].x_[2]=decstrx[tmp_tt].ptc.z;
                        tmp_ttx = tmp_ttx -1;
                    }else 
                    {
                        fin_mat[3*(tmr_numx-tmp_tttx)+1].x_[0]=decstrx[tmp_tt].x;
                        fin_mat[3*(tmr_numx-tmp_tttx)+1].x_[1]=decstrx[tmp_tt].y;
                        fin_mat[3*(tmr_numx-tmp_tttx)+1].x_[2]=decstrx[tmp_tt].z;
                        fin_mat[3*(tmr_numx-tmp_tttx)+0].x_[0]=decstrx[tmp_tt].ptn.x;
                        fin_mat[3*(tmr_numx-tmp_tttx)+0].x_[1]=decstrx[tmp_tt].ptn.y;
                        fin_mat[3*(tmr_numx-tmp_tttx)+0].x_[2]=decstrx[tmp_tt].ptn.z;
                        fin_mat[3*(tmr_numx-tmp_tttx)+2].x_[0]=decstrx[tmp_tt].ptc.x;
                        fin_mat[3*(tmr_numx-tmp_tttx)+2].x_[1]=decstrx[tmp_tt].ptc.y;
                        fin_mat[3*(tmr_numx-tmp_tttx)+2].x_[2]=decstrx[tmp_tt].ptc.z; 
                        tmp_tttx=tmp_tttx-1;                       
                    }
                    tmp_tt=tmp_tt+1;
                } 
 
            tmp_tt=0;
            for(int j=0;j<tmr_numx;j++)
            {
                fin_maty[3*(j+2)+0]=fin_mat[3*j+0];
                fin_maty[3*(j+2)+1]=fin_mat[3*j+1];
                fin_maty[3*(j+2)+2]=fin_mat[3*j+2];
                tmp_tt = tmp_tt+ 1;
            }
        //    cout<<"gggx"<<endl;
               
        } else
        {
            angle_rotate = (2.0*randf0and1()-1.0)*ang; // rotate angle
            GroupRotationpid(fin_maty[1].x_,fin_maty[fin_maty.size()-2].x_,angle_rotate,fin_maty,1,fin_maty.size()-2);            
        }


            dis_p=0.0;
         //   tp_clashx= false;
            new_vwd = 0.0;
            new_clash = 0.0;
            for(int jj=0;jj<fin_maty.size();jj++)
            {
                if(fin_maty[jj].elt_=="CA")
                {
        //            float VR1= GetVdwRadius(tmp_mat[jj].elt_);
        //            string VR1 = tmp_mat[jj].elt_;
                    for(int t=0;t<pointsBx.size();t++)
                    {
                        if(pointsBx[t].elt_=="CA")
                        {
            //                if((jj+rand_a)==t) continue;
                            if(abs(jj+3*rand_a-t)<1) continue;
            //                if(t>3*rand_a && t<(3*rand_a+tmp_mat.size())) continue;
            //                float VR2= GetVdwRadius(pointsBx[t].elt_);
            //                string VR2 = pointsBx[t].elt_ ;
                            dis_p = Distance_point(fin_maty[jj].x_,pointsBx[t].x_);
                            //if(dis_p < (VR1+VR2))
                            if(dis_p < 3.75)
                            {
                                new_clash = new_clash + 1.0/ sqrt(dis_p) ;
             //                   old_clash = old_clash + exp( 3.75-dis_p) ;
                            }                        
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
    //            cout<<"tp_clashx: " <<tp_clashx<<endl;
    //            if(tp_clashx) continue;  
    //            cout<<" 3 "<<endl;                 
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
                atm_idy=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
            bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
            bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;

            for(int j=0;j<fin_maty.size();j++)
            {
                if(fin_maty[j].elt_ == "CA")
                {
                    vector<float> cartX1;
                    vector<float> fracX1; 
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
                    atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);                 
                    for(int z=1;z<=theDensityMap.density.u3();z++)
                    {
                        atm_j[2] = z;
                        del_ij[2] =(atm_idy[tm_i1][2]-atm_j[2])/theDensityMap.grid[2];
                        if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                        if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                        del_ij[0] = del_ij[1] = 0.0;
                        vector<float> frac_tmpz;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpz);
                        if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;               
                        for(int y=1;y<=theDensityMap.density.u2();y++)
                        {
                            atm_j[1] = y;
                            del_ij[1] = (atm_idy[tm_i1][1] - atm_j[1])/theDensityMap.grid[1];
                            // wrap-around??
                            if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                            if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                            del_ij[0] = 0.0;
                            vector<float> frac_tmpy;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpy);
                            if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;                      
                            for(int x=1;x<=theDensityMap.density.u1();x++)
                            {
                                atm_j[0] = x;
                                del_ij[0] = (atm_idy[tm_i1][0] - atm_j[0])/theDensityMap.grid[0];
                                // wrap-around??
                                if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                vector<float> cart_del_ij2;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ij,cart_del_ij2);
                                float d2 = square_len(cart_del_ij2);
                                if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;
                            
                                float atm = C*exp(-k*d2);
                                float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                float inv_msk = 1.0/(1.0+sigmoid_msk);
                //                rhoC2(x,y,z) += atm;
                //                rhoC2(x,y,z) += atm;
                                rhoC2(x,y,z) += atm;
                                inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);


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
                    tm_i1 = tm_i1 + 1;
                }
            }
    //        rhoC2 = rhoC01;
            max_rhc = 0.0;
            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                if(rhoC2[x]>max_rhc) max_rhc = rhoC2[x] ;
            }     
    /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
        //        cout<<theDensityMap.density[x]<<" ";
                if(rhoC2[x]<(3.0/5.0)*max_rhc) rhoC2[x]=0; 
            } */
//            rhoC2 = rhoC01;                           
            sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
            sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
            clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
    /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
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
                }   */
            for(int x=bond0min;x<=bond0max;x++)
            {
               for(int y=bond1min;y<=bond1max;y++)
                {
                    for(int z=bond2min;z<=bond2max;z++)
                    {
           //            tmp_den=0.0;
            //            tmp_rhc =0.0;
            //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
            //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
            //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        clc_x2 = rhoC2(x,y,z);
            //            clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density(x,y,z);
            //            obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask2(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
            //            eps_x2 = 1.0;

                        // SMOOTHED
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
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }

                new_CC1 =1.0-CC_i2; 
        //        cout<< "new_CC: "<<new_CC<<" ";
            //    new_dE = 30*new_CC + new_vwd/((pnumx-1)*(rand_b-rand_a+1));
            //    new_dE = 0.8*new_CC + 0.5*new_clash/abs(rand_b-rand_a+1) + 0.5*new_Ehbond/abs(rand_b-rand_a+1);
                new_dE = 0.5*new_CC1 ;//+ 0.5*new_clash/abs(rand_b-rand_a+1);
            //    cout<<"new coor:  "<<tmp_mat[5].x_[0]<<" "<<tmp_mat[5].x_[1]<<" "<<tmp_mat[5].x_[2]<<endl;
                cout<<"old: "<< old_CC1<<" " <<old_clash<<endl;
                cout<<"new: "<< new_CC1<<" " <<new_clash<<endl;
        //        cout<<"old new dE: "<<old_dE<<" "<<new_dE<<" "<<dE<<endl;
                dE = new_dE - old_dE;
                if(new_dE <old_dE)
                {
                    int tmp_tmx=0;
                    for(int j=rand_a;j<=rand_b;j++)
                    {
                        pointsBx[3*j+0] = fin_maty[3*tmp_tmx+0];
                        pointsBx[3*j+1] = fin_maty[3*tmp_tmx+1];
                        pointsBx[3*j+2] = fin_maty[3*tmp_tmx+2];
                    //    bb[j]=bby[tmp_tmx];
                        tmp_tmx = tmp_tmx +1;
                    }
                //    rhoC0 =rhoC01;
                    cout<<" new_dE RRRRRRRRRRRR: "<<new_dE<<endl; 
                }                
                else
                {
                //    float tmpx=rand()/double(RAND_MAX);
                    float tmpx=randf0and1();
                    float mc_v = exp(-dE/(KT));
                    if(tmpx < mc_v)
                    {
                        int tmp_tmx=0;
                        for(int j=rand_a;j<=rand_b;j++)
                        {
                            pointsBx[3*j+0] = fin_maty[3*tmp_tmx+0];
                            pointsBx[3*j+1] = fin_maty[3*tmp_tmx+1];
                            pointsBx[3*j+2] = fin_maty[3*tmp_tmx+2];
                        //    bb[j]=bby[tmp_tmx];
                            tmp_tmx = tmp_tmx +1;
                        }
                    //    rhoC0 =rhoC01;
                        cout<<" XXX new_dE: "<<new_dE<<endl; 
                    }       
                }                         
    }

    if(dm_beg>3)
    {
        vector<int > beg_end;
        getLinLen(dm_beg,pointsBx,bb,beg_end);
        int Ni = beg_end[0]-1;
        int Ci = beg_end[1]+1;
        int LinLen = Ci - Ni + 1;  
        cout<<"The linker length is " <<LinLen <<" (from residue " << Ni<<" to "<<Ci<<")."<<endl; 

        int numseq=abs(Ci - Ni + 1);
//        if((Ci - Ni + 1)<3) continue;
        vector<point3f> decstr(LinLen);

            int tmp_tt=0;
            int tmp_ttx=0;
            for(int j=Ni;j<=Ci;j++)
            {
                decstr[tmp_tt].ss2=numtoss(bb[j].sst);
                decstr[tmp_tt].x=pointsBx[3*j+1].x_[0];
                decstr[tmp_tt].y=pointsBx[3*j+1].x_[1];
                decstr[tmp_tt].z=pointsBx[3*j+1].x_[2];
                decstr[tmp_tt].ptn.x=pointsBx[3*j+0].x_[0];
                decstr[tmp_tt].ptn.y=pointsBx[3*j+0].x_[1];
                decstr[tmp_tt].ptn.z=pointsBx[3*j+0].x_[2];
                decstr[tmp_tt].ptc.x=pointsBx[3*j+2].x_[0];
                decstr[tmp_tt].ptc.y=pointsBx[3*j+2].x_[1];
                decstr[tmp_tt].ptc.z=pointsBx[3*j+2].x_[2];
                tmp_tt=tmp_tt+1;
            }

        mcfragsweepLMP2(decstr,numseq,1,LinLen-2);
            tmp_tt=0;
            tmp_ttx=0;
            for(int j=Ni;j<=Ci;j++)
            {
                pointsBx[3*j+1].x_[0]=decstr[tmp_tt].x;
                pointsBx[3*j+1].x_[1]=decstr[tmp_tt].y;
                pointsBx[3*j+1].x_[2]=decstr[tmp_tt].z;
                pointsBx[3*j+0].x_[0]=decstr[tmp_tt].ptn.x;
                pointsBx[3*j+0].x_[1]=decstr[tmp_tt].ptn.y;
                pointsBx[3*j+0].x_[2]=decstr[tmp_tt].ptn.z;
                pointsBx[3*j+2].x_[0]=decstr[tmp_tt].ptc.x;
                pointsBx[3*j+2].x_[1]=decstr[tmp_tt].ptc.y;
                pointsBx[3*j+2].x_[2]=decstr[tmp_tt].ptc.z;
                tmp_tt=tmp_tt+1;
            }

    /*    if(LinLen >=3)
        {
            vector<poseCoord> fra_tmp;
            for(int tt=Ni;tt<=Ci;tt++)
            {
        //        fra_tmp.push_back(pointsBx[3*tt+0]);
                fra_tmp.push_back(pointsBx[3*tt+1]);
        //        fra_tmp.push_back(pointsBx[3*tt+2]);
            }
            connectgap(fra_tmp);
            int tmpee=0;
            for(int tt=Ni;tt<=Ci;tt++)
            {
                pointsBx[3*tt+0].x_[0]=pointsBx[3*tt+0].x_[0]+pointsBx[3*tt+1].x_[0]-fra_tmp[tmpee].x_[0];
                pointsBx[3*tt+0].x_[1]=pointsBx[3*tt+0].x_[1]+pointsBx[3*tt+1].x_[1]-fra_tmp[tmpee].x_[1];
                pointsBx[3*tt+0].x_[2]=pointsBx[3*tt+0].x_[2]+pointsBx[3*tt+1].x_[2]-fra_tmp[tmpee].x_[2];
                pointsBx[3*tt+1]=fra_tmp[tmpee];
                pointsBx[3*tt+2].x_[0]=pointsBx[3*tt+2].x_[0]+pointsBx[3*tt+1].x_[0]-fra_tmp[tmpee].x_[0];
                pointsBx[3*tt+2].x_[1]=pointsBx[3*tt+2].x_[1]+pointsBx[3*tt+1].x_[1]-fra_tmp[tmpee].x_[1];
                pointsBx[3*tt+2].x_[2]=pointsBx[3*tt+2].x_[2]+pointsBx[3*tt+1].x_[2]-fra_tmp[tmpee].x_[2];
                tmpee = tmpee +1;
            }            
        }   */

    }
    if(dm_end<(pnum-4) && dm_end>4)
    {
        vector<int > beg_end;
        getLinLen(dm_end,pointsBx,bb,beg_end);
        int Ni = beg_end[0]-1;
        int Ci = beg_end[1]+1;
        int LinLen = Ci - Ni + 1;  
        cout<<"The linker length is " <<LinLen <<" (from residue " << Ni<<" to "<<Ci<<")."<<endl; 

        vector<point3f> decstrx(LinLen);
        int numseq=abs(Ci - Ni + 1);
        int tmp_tt=0;
    //    if((rand_b-rand_a-3)<0) continue;
            for(int j=Ni;j<=Ci;j++)
            {
                decstrx[tmp_tt].ss2=numtoss(bb[j].sst);
                decstrx[tmp_tt].x=pointsBx[3*j+1].x_[0];
                decstrx[tmp_tt].y=pointsBx[3*j+1].x_[1];
                decstrx[tmp_tt].z=pointsBx[3*j+1].x_[2];
                decstrx[tmp_tt].ptn.x=pointsBx[3*j+0].x_[0];
                decstrx[tmp_tt].ptn.y=pointsBx[3*j+0].x_[1];
                decstrx[tmp_tt].ptn.z=pointsBx[3*j+0].x_[2];
                decstrx[tmp_tt].ptc.x=pointsBx[3*j+2].x_[0];
                decstrx[tmp_tt].ptc.y=pointsBx[3*j+2].x_[1];
                decstrx[tmp_tt].ptc.z=pointsBx[3*j+2].x_[2];
                tmp_tt=tmp_tt+1;
            }            
        mcfragsweepLMP2(decstrx,numseq,1,LinLen-2);
        tmp_tt=0;
        for(int j=Ni;j<=Ci;j++)
        {
            pointsBx[3*j+1].x_[0]=decstrx[tmp_tt].x;
            pointsBx[3*j+1].x_[1]=decstrx[tmp_tt].y;
            pointsBx[3*j+1].x_[2]=decstrx[tmp_tt].z;
            pointsBx[3*j+0].x_[0]=decstrx[tmp_tt].ptn.x;
            pointsBx[3*j+0].x_[1]=decstrx[tmp_tt].ptn.y;
            pointsBx[3*j+0].x_[2]=decstrx[tmp_tt].ptn.z;
            pointsBx[3*j+2].x_[0]=decstrx[tmp_tt].ptc.x;
            pointsBx[3*j+2].x_[1]=decstrx[tmp_tt].ptc.y;
            pointsBx[3*j+2].x_[2]=decstrx[tmp_tt].ptc.z;
            tmp_tt=tmp_tt+1;
        } 

   /*     if(LinLen >=3)
        {
            vector<poseCoord> fra_tmp;
            for(int tt=Ni;tt<=Ci;tt++)
            {
    //            fra_tmp.push_back(pointsBx[3*tt+0]);
                fra_tmp.push_back(pointsBx[3*tt+1]);
    //            fra_tmp.push_back(pointsBx[3*tt+2]);
            }
            connectgap(fra_tmp);
            int tmpee=0;
            for(int tt=Ni;tt<=Ci;tt++)
            {
                pointsBx[3*tt+0].x_[0]=pointsBx[3*tt+0].x_[0]+pointsBx[3*tt+1].x_[0]-fra_tmp[tmpee].x_[0];
                pointsBx[3*tt+0].x_[1]=pointsBx[3*tt+0].x_[1]+pointsBx[3*tt+1].x_[1]-fra_tmp[tmpee].x_[1];
                pointsBx[3*tt+0].x_[2]=pointsBx[3*tt+0].x_[2]+pointsBx[3*tt+1].x_[2]-fra_tmp[tmpee].x_[2];
                pointsBx[3*tt+1]=fra_tmp[tmpee];
                pointsBx[3*tt+2].x_[0]=pointsBx[3*tt+2].x_[0]+pointsBx[3*tt+1].x_[0]-fra_tmp[tmpee].x_[0];
                pointsBx[3*tt+2].x_[1]=pointsBx[3*tt+2].x_[1]+pointsBx[3*tt+1].x_[1]-fra_tmp[tmpee].x_[1];
                pointsBx[3*tt+2].x_[2]=pointsBx[3*tt+2].x_[2]+pointsBx[3*tt+1].x_[2]-fra_tmp[tmpee].x_[2];
                tmpee = tmpee +1;
            }            
        }    */         
    }     

    cout<< "TTTTTT"<<endl;
    for(int jjj=-1;jjj>=0;jjj--)
    {
    vector<poseCoord > tmp_mat;
    vector<float> cartX2;
    vector<float> fracX2;
    float sumis=0.0; // intersection
    float sumu=0.0; // union
//    string elt_i;
    int allnum=500;
    int pnum_tmp=0;
    int dm_s=0; // start point for each domain
    int dm_e=0; // end point for each domain
    if(jjj==0) 
    {
        allnum= (randop-1-0+1)*allnum;
        pnum_tmp=randop-1-0+1;
        dm_s=0;
        dm_e= randop;
    }
    if(jjj==1)
    {
        allnum = (pnum-randop)*allnum;
        pnum_tmp= pnum-randop;
        dm_s=87;
        dm_e=pnum;
    }
//    int allnum=1;
    for(int tt=0;tt<allnum;tt++)
    {       
//        int randpx = 36;
        float KT0 =0.001;
        float KTn =0.0001;
        float KTx= pow(float(KTn/KT0),float(float(tt)/float(5)));
//        KT = KT0*float(KTx); 
        KT=0.001;
//        float angle0 =180.0;
        float angle0 =3.0;
//        float anglen =10;
        float anglen =1.0;
        float anglex= pow(float(anglen/angle0),float(float(tt)/float(5)));
//        ang = angle0*float(anglex);
        ang=2.0;
//        int randp = rand()%15+30;
        int recT = 0;   
        int rec_t = 0 ;  
        int rand_a = 0;
        int rand_b = 0; 
        int flagx=0 ;

            vector<poseCoord>().swap(tmp_mat);
            //    int randtx = rand()%pnum_tmp+dm_s;
                int randtx = randIntCustom(0,pnum_tmp-1)+dm_s;
                int randt=0;
                if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                {
                //    randt=rand()%4+6;
                    randt=randIntCustom(0,3)+6;
                    int jmx=randtx,jnx=randtx;
                    for(int j=0;j<randt;j++)
                    {
                        jnx=randtx+j;
                        if(jnx>=(dm_e-1)) break;
                        if(numtoss(bb[jnx].sst)=='C') break;
                    }
                    if((jnx-randtx+1)<randt)
                    {
                        for(int j=0;j<(randt-(jnx-randtx));j++)
                        {
                            jmx=randtx-j;
                            if(jmx<=dm_s) break;
                            if(numtoss(bb[jmx].sst)=='C') break;
                        }
                    }
                    rand_a=jmx;
                    rand_b=jnx;
                } else
                {
                //    randt = rand()%4+3;
                    randt = randIntCustom(0,3)+3;
                    int jmx=randtx,jnx=randtx;
                    for(int j=0;j<randt;j++)
                    {
                        jnx=randtx+j;
                        if(jnx>=(dm_e-1)) break;
                        if(numtoss(bb[jnx].sst)!='C') break;
                    }
                    if((jnx-randtx+1)<randt)
                    {
                        for(int j=0;j<(randt-(jnx-randtx));j++)
                        {
                            jmx=randtx-j;
                            if(jmx<=dm_s) break;
                            if(numtoss(bb[jmx].sst)!='C') break;
                        }
                    }
                    rand_a=jmx;
                    rand_b=jnx;                    
                }
            old_dE = 0.0;
            if(rand_b>(dm_e-1)||rand_a<dm_s) continue;
            int randtx_x=randtx;
            if(randtx>=(dm_e-4)||randtx<=dm_s+3)
            {                  
                if(randtx>=(dm_e-4))
                {
                    if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                    {
                        int jmx=0,jnx=0;
                        if(randtx>=(dm_e-4))
                        {
                            for(int j=0;j<pnum_tmp;j++)
                            {
                                jmx=randtx-j;
                                if(jmx<=dm_s)
                                {
                                    jmx=dm_s;
                                    break;
                                }
                                if(numtoss(bb[jmx].sst)=='C') break;
                            }
                            randtx_x = jmx;              
                        }                   
                    }  
                    rand_a=randtx_x;
                    rand_b=dm_e-1;
                    for(int j=rand_a;j<=rand_b;j++)
                    {
                        tmp_mat.push_back(pointsBx[3*j+0]);
                        tmp_mat.push_back(pointsBx[3*j+1]);
                        tmp_mat.push_back(pointsBx[3*j+2]);
                    }
                } 
                if(randtx<=dm_s+3)
                {
                    if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                    {
                        int jmx=0,jnx=0;
                        if(randtx<=dm_s+3)
                        {
                            for(int j=0;j<pnum_tmp;j++)
                            {
                                jnx=randtx+j;
                                if(jnx>=dm_e-1)
                                {
                                    jnx=dm_e-1;
                                    break;
                                }
                                if(numtoss(bb[jnx].sst)=='C') break;
                            }
                            randtx_x = jnx;            
                        }                   
                    }                     
                    rand_a=0+dm_s;
                    rand_b=randtx_x;
                    for(int j=rand_a;j<=rand_b;j++)
                    {
                        tmp_mat.push_back(pointsBx[3*j+0]);
                        tmp_mat.push_back(pointsBx[3*j+1]);
                        tmp_mat.push_back(pointsBx[3*j+2]);
                    }                    
                }               
            } else
            {
                if(rand_a>rand_b)
                {
                    int tmp_rand = rand_b;
                    rand_b = rand_a;
                    rand_a = tmp_rand;
                }                
//                if( abs(rand_a -rand_b+1)<4) continue;
                for(int j=rand_a;j<=rand_b;j++)
                {
                    tmp_mat.push_back(pointsBx[3*j+0]);
                    tmp_mat.push_back(pointsBx[3*j+1]);
                    tmp_mat.push_back(pointsBx[3*j+2]);
                }
            }
            

            float old_Ehbond=0.0;
            int numseq=abs(rand_b-rand_a+1);
            if(numseq<3) continue;
            vector<point3f> decstr(numseq);

            if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
            {
                tmp_tt=0;
                for(int j=rand_a;j<=rand_b;j++)
                {
                    decstr[tmp_tt].ss2=numtoss(bb[j].sst);
                    decstr[tmp_tt].x=tmp_mat[3*tmp_tt+1].x_[0];
                    decstr[tmp_tt].y=tmp_mat[3*tmp_tt+1].x_[1];
                    decstr[tmp_tt].z=tmp_mat[3*tmp_tt+1].x_[2];
                    decstr[tmp_tt].ptn.x=tmp_mat[3*tmp_tt+0].x_[0];
                    decstr[tmp_tt].ptn.y=tmp_mat[3*tmp_tt+0].x_[1];
                    decstr[tmp_tt].ptn.z=tmp_mat[3*tmp_tt+0].x_[2];
                    decstr[tmp_tt].ptc.x=tmp_mat[3*tmp_tt+2].x_[0];
                    decstr[tmp_tt].ptc.y=tmp_mat[3*tmp_tt+2].x_[1];
                    decstr[tmp_tt].ptc.z=tmp_mat[3*tmp_tt+2].x_[2];
                    tmp_tt=tmp_tt+1;
                } 
                old_Ehbond=energyhbondcanc(decstr,numseq);
            }
        //    cout<<"old_Ehbond: "<<old_Ehbond<<" ";

            // ss information
    /*        float old_ss_E=0.0;
            int CA_num=0;
            vector<vector<float> > old_CA_cor;
            int tm_i1=0;               
            for(int j=0;j<tmp_mat.size();j++)
            { 
                if(tmp_mat[j].elt_ == "CA")
                {
                    old_CA_cor.push_back(tmp_mat[j].x_);
                    CA_num = CA_num+1;
                    tm_i1=tm_i1+1; 
                }                  
            }
            vector<string > old_ss_inf(CA_num);
            make_secx(old_CA_cor,CA_num,old_ss_inf);
            int ss_ab_num=0; // helix and b strand number
            int ss_c_num =0; // coil number
//                float old_ss_ab=0.0;
//                float old_ss_c=0.0;                
            for(int j=0;j<old_ss_inf.size();j++)
            {
                if(old_ss_inf[j]=='H')
                {
                    old_ss_E = old_ss_E +0.1;
//                        old_ss_ab = old_ss_ab + 0.4;
//                        ss_ab_num = ss_ab_num +1 ;
                }
                else if(old_ss_inf[j]=='E')
                {
                    old_ss_E = old_ss_E +0.1;
//                        old_ss_ab = old_ss_ab + 0.4;
//                        ss_ab_num = ss_ab_num +1;
                }
                else{
                    old_ss_E = old_ss_E +0.2;
//                        old_ss_c = old_ss_c + 0.5;
//                        ss_c_num = ss_c_num +1;
                }                    
            }
//                old_ss_E = (old_ss_ab*ss_ab_num+old_ss_c*ss_c_num)/(ss_ab_num+ss_c_num);
            old_ss_E = old_ss_E/(CA_num); */

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
            int tm_i1=0;          
                del_ij=vector<float>(3,0.0);
                atm_j=vector<float>(3,0.0);
                tm_i1 =0;
                vector<vector<float>> atm_idy((rand_b-rand_a+1),vector<float>(3,0.0));
            bond0max=0.0;bond1max=0.0;bond2max=0.0;
            bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;                
    //            atm_idy=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
            for(int j=0;j<tmp_mat.size();j++)
            {
                if(tmp_mat[j].elt_ == "CA")
                {
                    vector<float> cartX1;
                    vector<float> fracX1; 
                //    elt_i = tmp_mat[3*tm_i1+1].elt_;
                    elt_i = tmp_mat[j].elt_;
                    elt_i = elt_i[0];
                    OneGaussianScattering sig_j = get_A( elt_i );
                    k = sig_j.k( theDensityMap.effectiveB );
                    C = sig_j.C( k );
                    if ( C < 1e-6 ) continue; 

            //        cartX1 = tmp_mat[3*tm_i1+1].x_;
                    cartX1 = tmp_mat[j].x_;
            //        tm_i1 = tm_i1 + 1;
                    MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                    atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);                 
                    for(int z=1;z<=theDensityMap.density.u3();z++)
                    {
                        atm_j[2] = z;
                        del_ij[2] =(atm_idy[tm_i1][2]-atm_j[2])/theDensityMap.grid[2];
                        if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                        if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                        del_ij[0] = del_ij[1] = 0.0;
                        vector<float> frac_tmpz;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpz);
                        if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;               
                        for(int y=1;y<=theDensityMap.density.u2();y++)
                        {
                            atm_j[1] = y;
                            del_ij[1] = (atm_idy[tm_i1][1] - atm_j[1])/theDensityMap.grid[1];
                            // wrap-around??
                            if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                            if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                            del_ij[0] = 0.0;
                            vector<float> frac_tmpy;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpy);
                            if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;                      
                            for(int x=1;x<=theDensityMap.density.u1();x++)
                            {
                                atm_j[0] = x;
                                del_ij[0] = (atm_idy[tm_i1][0] - atm_j[0])/theDensityMap.grid[0];
                                // wrap-around??
                                if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                vector<float> cart_del_ij2;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ij,cart_del_ij2);
                                float d2 = square_len(cart_del_ij2);
                                if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;
                            
                                float atm = C*exp(-k*d2);
                                float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                float inv_msk = 1.0/(1.0+sigmoid_msk);
                //                rhoC2(x,y,z) += atm;
                //                rhoC2(x,y,z) += atm;
                                rhoC2(x,y,z) += atm;
                 //               if(rhoC01(x,y,z)<0) rhoC01(x,y,z)=0;
                                inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);

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
                    tm_i1 = tm_i1 + 1;
                }
            }

 //           for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
 //               rhoC2[x] = rhoC0[x];
//            }            
            float max_rhc = 0.0;
            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                if(rhoC2[x]>max_rhc) max_rhc = rhoC2[x] ;
            }
    /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
        //        cout<<theDensityMap.density[x]<<" ";
                if(rhoC2[x]<(3.0/5.0)*max_rhc) rhoC2[x]=0; 
            }    */
//            rhoC2 = rhoC0;         
            float sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
            float sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
            float clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
/*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
//                tmp_den = 0.0; 
//                tmp_rhc=0.0;
  //              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
//                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                // fetch this point
                clc_x2 = rhoC2[x];
//                clc_x2 = tmp_rhc;
                obs_x2 = theDensityMap.density[x];
//                obs_x2 = tmp_den;
                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
            //    eps_x2 = 1.0;

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
            }  */
            for(int x=int(bond0min);x<=bond0max;x++)
            {
               for(int y=int(bond1min);y<=bond1max;y++)
                {
                    for(int z=int(bond2min);z<=bond2max;z++)
                    {
        //                tmp_den = 0.0; 
        //                tmp_rhc=0.0;
          //              if(theDensityMap.density[x]>(2.0/3.0)*max_den) tmp_den=theDensityMap.density[x];
        //                if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
        //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        clc_x2 = rhoC2(x,y,z);
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density(x,y,z);
        //                obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask2(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                    //    eps_x2 = 1.0;

                        // SMOOTHED
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
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }           
    /*        sumis=0.0;
            sumu=0.0;
            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                tmp_den = 0.0; 
                tmp_rhc=0.0;
                if(theDensityMap.density[x]>(3.0/5.0)*max_den) tmp_den=1.0;
                if(rhoC2[x]>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                // fetch this point
                sumis=sumis+tmp_den*tmp_rhc;
                sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
            }
            if ( sumu == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = sumis/sumu;
            }            */
//            vector<vector<float> > bondx(3,vector<float>(2,0.0));
   /*         float bond0max=0.0,bond1max=0.0,bond2max=0.0;
            float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;
//            vector<vector<float>> atm_idy((rand_b-rand_a+1),vector<float>(3,0.0));
            atm_idy=vector<vector<float>>((rand_b-rand_a+1),vector<float>(3,0.0));
            tm_i1=0; 
            for(int j=0;j<tmp_mat.size();j++)
            {
 //               tmp_mat.push_back(pointsBx[j]);
                if(tmp_mat[j].elt_ == "CA")
                {
                    vector<float> cartX1;
                    vector<float> fracX1;
    //                cartX1 = tmp_mat[3*tm_i1+1].x_;
                    cartX1 = tmp_mat[j].x_;
                    MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                    atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);
                    if(atm_idy[tm_i1][0]>bond0max) bond0max = atm_idy[tm_i1][0];
                    if(atm_idy[tm_i1][1]>bond1max) bond1max = atm_idy[tm_i1][1];
                    if(atm_idy[tm_i1][2]>bond2max) bond2max = atm_idy[tm_i1][2];
                    if(atm_idy[tm_i1][0]<bond0min) bond0min = atm_idy[tm_i1][0];
                    if(atm_idy[tm_i1][1]<bond1min) bond1min = atm_idy[tm_i1][1];
                    if(atm_idy[tm_i1][2]<bond2min) bond2min = atm_idy[tm_i1][2];
                    tm_i1 = tm_i1 + 1;
                }
            }  */
    /*        sumis=0.0;
            sumu=0.0;            
            for(int x=int(bond0min);x<=bond0max;x++)
            {
               for(int y=int(bond1min);y<=bond1max;y++)
                {
                    for(int z=int(bond2min);z<=bond2max;z++)
                    {
                        tmp_den = 0.0; 
                        tmp_rhc=0.0;
                        if(theDensityMap.density(x,y,z)>=(3.0/5.0)*max_den) tmp_den=1.0;
                        if(rhoC2(x,y,z)>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
        //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        sumis=sumis+tmp_den*tmp_rhc;
                        sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
                    }
                }
            }          
            if ( sumu == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = sumis/sumu;
            }  */

/*            float bond0max=0.0,bond1max=0.0,bond2max=0.0;
            float bond0min=1000.0,bond1min=1000.0,bond2min=1000.0;
            vector<vector<float>> atm_idy((rand_b-rand_a+1),vector<float>(3,0.0));
            tm_i1=0; 
            for(int j=0;j<tmp_mat.size();j++)
            {
 //               tmp_mat.push_back(pointsBx[j]);
                if(tmp_mat[j].elt_ == "CA")
                {
                    vector<float> cartX1;
                    vector<float> fracX1;
    //                cartX1 = tmp_mat[3*tm_i1+1].x_;
                    cartX1 = tmp_mat[j].x_;
                    MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                    atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);
                    if(atm_idy[tm_i1][0]>bond0max) bond0max = atm_idy[tm_i1][0];
                    if(atm_idy[tm_i1][1]>bond1max) bond1max = atm_idy[tm_i1][1];
                    if(atm_idy[tm_i1][2]>bond2max) bond2max = atm_idy[tm_i1][2];
                    if(atm_idy[tm_i1][0]<bond0min) bond0min = atm_idy[tm_i1][0];
                    if(atm_idy[tm_i1][1]<bond1min) bond1min = atm_idy[tm_i1][1];
                    if(atm_idy[tm_i1][2]<bond2min) bond2min = atm_idy[tm_i1][2];
                    tm_i1 = tm_i1 + 1;
                }
            } 
            sumis=0.0;
            sumu=0.0;            
            for(int x=int(bond0min);x<=bond0max;x++)
            {
               for(int y=int(bond1min);y<=bond1max;y++)
                {
                    for(int z=int(bond2min);z<=bond2max;z++)
                    {
                        clc_x2 = rhoC2(x,y,z);
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density(x,y,z);
        //                obs_x2 = tmp_den;
        //                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                        eps_x2 = 1.0;

                        // SMOOTHED
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
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }  */

/*                for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                    clc_x2 = rhoC2[x];
                    obs_x2 = theDensityMap.density[x];
                    eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal

                    // SMOOTHED
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                }
                sumO_i2 = sumO_i2/vol_i2;
                sumC_i2 = sumC_i2/vol_i2;
                sumCO_i2 = sumCO_i2/vol_i2;
                sumC2_i2 = sumC2_i2/vol_i2;
                sumO2_i2 = sumO2_i2/vol_i2;
                varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2);
                varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2) ;
                if ( varC_i2 == 0 || varO_i2 == 0 ) {
                    CC_i2 = 0;
                } else {
                    CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2) / sqrt( varC_i2 * varO_i2 );
                } */
        /*        for(int x=bond0min;x<=bond0max;x++)
                {
                   for(int y=bond1min;y<=bond1max;y++)
                    {
                        for(int z=bond2min;z<=bond2max;z++)
                        {
                            clc_x2 = rhoC2(x,y,z);
            //                clc_x2 = tmp_rhc;
                            obs_x2 = theDensityMap.density(x,y,z);
            //                obs_x2 = tmp_den;
            //                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                            eps_x2 = 1.0;

                            // SMOOTHED
                            vol_i2 += eps_x2;
                            sumC_i2 += (eps_x2*clc_x2*clc_x2);
                            sumO_i2 += (eps_x2*obs_x2*obs_x2);
                            sumCO_i2 += (eps_x2*clc_x2*obs_x2);
                        }
                    }
                }          
                if ( sumC_i2 == 0 || sumO_i2 == 0 ) {
                    CC_i2 = 0;
                } else {
                    CC_i2 = sumCO_i2/(sqrt(sumC_i2+sumO_i2)); ;
                }             */
    /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x )
                {
                    eps_x2 = 1-inv_rho_mask2[x];
                    vol_i2   += eps_x2;
                    sumC_i2 += (eps_x2*rhoC2[x]*rhoC2[x]);
                    sumO_i2 += (eps_x2*theDensityMap.density[x]*theDensityMap.density[x]);
                    sumCO_i2 += (eps_x2*rhoC2[x]*theDensityMap.density[x]);
                }
        //        sumCO_i2 = sumCO_i2/vol_i2;
        //        sumO_i2 = sumO_i2/vol_i2;
        //        sumC_i2 = sumC_i2/vol_i2;
                CC_i2 = sumCO_i2/(sqrt(sumC_i2*sumO_i2));  */
                old_CC =1.0- CC_i2; 
        //        cout<< "old_CC: "<<old_CC<<" ";

                float dis_p=0.0;
                bool tp_clashx= false;
                old_vwd =0.0;
                old_clash = 0.0;
                for(int jj=0;jj<tmp_mat.size();jj++)
                {
                    if(tmp_mat[jj].elt_=="CA")
                    {
//                        float VR1= GetVdwRadius(tmp_mat[jj].elt_);
                    //    string VR1 = tmp_mat[jj].elt_ ;
                    //    for(int t=0;t<pointsBx.size();t++)
                        for(int t=3*dm_s+0;t<3*dm_e;t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                //                if((jj+rand_a)==t) continue;
                                if(abs(jj+3*rand_a-t)<1) continue;
                       //         if(t>3*rand_a && t<(3*rand_a+tmp_mat.size())) continue;
//                                float VR2= GetVdwRadius(pointsBx[t].elt_);
                //                string VR2 = pointsBx[t].elt_ ;
                                dis_p = Distance_point(tmp_mat[jj].x_,pointsBx[t].x_);
//                                if(dis_p < (VR1+VR2))
                                if(dis_p < 3.75)
                                {
                                    old_clash = old_clash + 1.0/ sqrt(dis_p) ;
                //                    old_clash = old_clash + exp( 3.75-dis_p) ;
                                }
                //                old_vwd=old_vwd+GetVdwEg(VR1,VR2,dis_p);
                //                old_vwd=old_vwd+GetVdwEgC(VR1,VR2,dis_p);
            //                    cout<<"jj,t: "<<jj+rand_a<<" "<<t<<" dis: "<<VR1+VR2<<" "<<dis_p<<endl;
            //                    cout<<"vwdeg: "<<GetVdwEg(VR1,VR2,dis_p)<<endl;;
            //                    if(((dis_p)<2.0))
            //                    {
            //                        tp_clashx = true;
            //                        break;
            //                    } 
            //                    if(tp_clashx)
            //                    {
            //                        cout<<"jj,t: "<<tmp_mat[jj].elt_<<" "<<pointsBx[t].elt_<<" dis: "<<VR1+VR2<<" "<<dis_p<<endl;
            //                        tp_clashx=false;
            //                    }
            //                    cout<< "VR1+VR2,dis_p: "<<VR1+VR2<<" "<<dis_p<<endl;
                            }
                        }
                    }
    //                if(tp_clashx) break;
                }  

    //            cout<<"tp_clashx: " <<tp_clashx<<endl;
    //            if(tp_clashx) continue;  
    //            old_dE = old_CC + old_vwd;
    //            old_dE = 30*old_CC + old_vwd/((pnum-1)*(rand_b-rand_a+1));
           //     old_dE = 0.8*old_CC + 0.5*old_clash/abs(rand_b-rand_a+1)+0.5*old_Ehbond/abs(rand_b-rand_a+1);
                old_dE = 0.5*old_CC + 0.5*old_clash/abs(rand_b-rand_a+1)+0.5*old_Ehbond/abs(rand_b-rand_a+1);
           //     cout<< "old_dE: "<<old_dE<<" ";
           //     cout<<"old size: "<<tmp_mat.size()<<endl;
           //     cout<<"old coor:  "<<tmp_mat[5].x_[0]<<" "<<tmp_mat[5].x_[1]<<" "<<tmp_mat[5].x_[2]<<endl;



    //            rhoC01 = rhoC0;
                 
                
                // Rotate angle
        //        ang=180.0;      
//                vector<vector<float> > tmm;
//                for(int tc=0;tc<tmp_mat.size();tc++)
//                {
//                    tmm.push_back(tmp_mat[tc].x_);
//                }
//                GroupRotation(tmm[0],tmm[tmm.size()-1],angle_rotate,tmm);
//                for(int tc=0;tc<tmp_mat.size();tc++)
//                {
//                    tmp_mat[tc].x_ = tmm[tc];
//                }                
        //        GroupRotationp(tmp_mat[0].x_,tmp_mat[tmp_mat.size()-1].x_,angle_rotate,tmp_mat);
                if(randtx>=(dm_e-4)||randtx<=dm_s+3)
                {  // fragment in the N terminal or C terminal
            /*        if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                    {
                        int jmx=0,jnx=0;
                        if(randtx>(pnum-4))
                        {
                            for(int j=0;j<pnum;j++)
                            {
                                jmx=randtx-j;
                                if(jmx<=0)
                                {
                                    jmx=0;
                                    break;
                                }
                                if(numtoss(bb[jmx].sst)=='C') break;

                            //    if(numtoss(bb[jmx].sst)!='H'&&numtoss(bb[jmx].sst)!='E') break;
                            }
                            randtx = jmx;              
                        }
                        if(randtx<3)
                        {
                            for(int j=0;j<pnum;j++)
                            {
                                jnx=randtx+j;
                                if(jnx>=pnum-1)
                                {
                                    jnx=pnum-1;
                                    break;
                                }
                                if(numtoss(bb[jnx].sst)=='C') break;
                            //    if(numtoss(bb[jmx].sst)!='H'&&numtoss(bb[jmx].sst)!='E') break;
                            }
                            randtx = jnx;            
                        }                   
                    }  */
                    // Rotation axis 
                    float asin_theta=2.0*randf0and1()-1.0;
                    float acos_theta=sqrt(1.0-asin_theta*asin_theta);
                    float apha=2.0*PI*randf0and1();
                    float awx=acos_theta*cos(apha);
                    float awy=acos_theta*sin(apha);
                    float awz=asin_theta;
                    // Translation Vector
                    float t0=0.1;
                    float t1=(randf0and1()*2.0-1.0)*t0+0.0;
                    float t2=(randf0and1()*2.0-1.0)*t0+0.0;
                    float t3=(randf0and1()*2.0-1.0)*t0+0.0;

                    // Rotation matrix
                    float anggg=ang;
                    float angle_rotategg=(2.0*randf0and1()-1.0)*anggg; // rotate angle
                    float asin=sin(angle_rotategg);
                    float acos=cos(angle_rotategg);
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
                    int rand_point = randtx_x;
                    axyz[0]=pointsBx[3*rand_point+1].x_[0]; // CA atom
                    axyz[1]=pointsBx[3*rand_point+1].x_[1];
                    axyz[2]=pointsBx[3*rand_point+1].x_[2];     

                    int tmp3=0;
                    vector<poseCoord > fin_mat;
                    fin_mat = tmp_mat;
                    if(randtx>=(dm_e-4))
                    {
                        for(int j=1;j<tmp_mat.size();j++) // not first N
                        {
                            fin_mat[tmp3].x_[0]=t1+axyz[0]+(tmp_mat[j].x_[0]-axyz[0])*u[0][0]+(tmp_mat[j].x_[1]-axyz[1])*u[0][1]+(tmp_mat[j].x_[2]-axyz[2])*u[0][2];
                            fin_mat[tmp3].x_[1]=t2+axyz[1]+(tmp_mat[j].x_[0]-axyz[0])*u[1][0]+(tmp_mat[j].x_[1]-axyz[1])*u[1][1]+(tmp_mat[j].x_[2]-axyz[2])*u[1][2];
                            fin_mat[tmp3].x_[2]=t3+axyz[2]+(tmp_mat[j].x_[0]-axyz[0])*u[2][0]+(tmp_mat[j].x_[1]-axyz[1])*u[2][1]+(tmp_mat[j].x_[2]-axyz[2])*u[2][2];
                            tmp3=tmp3+1;
                        }
                    } else
                    {
                        for(int j=0;j<tmp_mat.size()-1;j++) // not last N
                        {
                            fin_mat[tmp3].x_[0]=t1+axyz[0]+(tmp_mat[j].x_[0]-axyz[0])*u[0][0]+(tmp_mat[j].x_[1]-axyz[1])*u[0][1]+(tmp_mat[j].x_[2]-axyz[2])*u[0][2];
                            fin_mat[tmp3].x_[1]=t2+axyz[1]+(tmp_mat[j].x_[0]-axyz[0])*u[1][0]+(tmp_mat[j].x_[1]-axyz[1])*u[1][1]+(tmp_mat[j].x_[2]-axyz[2])*u[1][2];
                            fin_mat[tmp3].x_[2]=t3+axyz[2]+(tmp_mat[j].x_[0]-axyz[0])*u[2][0]+(tmp_mat[j].x_[1]-axyz[1])*u[2][1]+(tmp_mat[j].x_[2]-axyz[2])*u[2][2];
                            tmp3=tmp3+1;
                        }                        
                    }
                    vector<poseCoord>().swap(tmp_mat);
                    tmp_mat=fin_mat;

                } else 
                {
                    angle_rotate = (2.0*randf0and1()-1.0)*ang; // rotate angle 
                    GroupRotationpid(tmp_mat[1].x_,tmp_mat[tmp_mat.size()-2].x_,angle_rotate,tmp_mat,1,tmp_mat.size()-2);
                }
//                cout<<" 1 ";
/*                int randxx = rand()%(rand_b-rand_a+1)+rand_a;
                vector<float> tmpcc(3,0.0);
                int randyy=randxx-rand_a;  // randxx and randyy is the same point
                float deta0 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.1+0.0;
                float deta1 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.1+0.0;
                float deta2 = (RandomDoubleX(0.0,1.0)*2.0-1.0)*0.1+0.0;
                tmpcc[0] = tmp_mat[randyy].x_[0] + deta0;
                tmpcc[1] = tmp_mat[randyy].x_[1] + deta1;
                tmpcc[2] = tmp_mat[randyy].x_[2] + deta2;
                float dis_s0=0.0,dis_s1=0.0;
                if(randxx==0) {
                    dis_s0 = Distance_point(tmpcc,pointsBx[randxx+1].x_);
                }
                else if(randxx==(pnum-1))
                {
                    dis_s0 = Distance_point(tmpcc,pointsBx[randxx-1].x_);
                }
                else {
                    if(randyy-1>=0) dis_s0 = Distance_point(tmpcc,tmp_mat[randyy-1].x_);
                    if(randyy+1 < tmp_mat.size()) dis_s0 = Distance_point(tmpcc,tmp_mat[randyy+1].x_);
                    if(randyy==(tmp_mat.size()-1)&& (randxx+1)<=(pnumx-1)) dis_s0 = Distance_point(tmpcc,pointsBx[randxx+1].x_);
                    if(randyy==0&&(randxx-1)>=0) dis_s0 = Distance_point(tmpcc,pointsBx[randxx-1].x_);
                } */
        //      if(dis_s0>3.70 && dis_s0<3.85 ) tmp_mat[randyy].x_ = tmpcc;   
        //        if(dis_s0>1.43001 && dis_s0<1.49999 && dis_s1<(1.33999) && dis_s1>(1.31001)) tmp_mat[randyy] = tmpcc;
//              cout<<" 2 ";

                dis_p=0.0;
                tp_clashx= false;
                new_vwd = 0.0;
                new_clash = 0.0;
                for(int jj=0;jj<tmp_mat.size();jj++)
                {
                    if(tmp_mat[jj].elt_=="CA")
                    {
            //            float VR1= GetVdwRadius(tmp_mat[jj].elt_);
            //            string VR1 = tmp_mat[jj].elt_;
                    //    for(int t=0;t<pointsBx.size();t++)
                        for(int t=3*dm_s+0;t<3*dm_e;t++)
                        {
                            if(pointsBx[t].elt_=="CA")
                            {
                //                if((jj+rand_a)==t) continue;
                                if(abs(jj+3*rand_a-t)<1) continue;
                //                if(t>3*rand_a && t<(3*rand_a+tmp_mat.size())) continue;
                //                float VR2= GetVdwRadius(pointsBx[t].elt_);
                //                string VR2 = pointsBx[t].elt_ ;
                                dis_p = Distance_point(tmp_mat[jj].x_,pointsBx[t].x_);
                                //if(dis_p < (VR1+VR2))
                                if(dis_p < 3.75)
                                {
                                    new_clash = new_clash + 1.0/ sqrt(dis_p) ;
                 //                   old_clash = old_clash + exp( 3.75-dis_p) ;
                                }                        
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
    //            cout<<"tp_clashx: " <<tp_clashx<<endl;
    //            if(tp_clashx) continue;  
    //            cout<<" 3 "<<endl;                 

                for ( int t=0; t<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++t ) {
                    rhoC2[t]=0.0;
                    inv_rho_mask2[t]=1.0;
                }
                del_ij=vector<float>(3,0.0);
                atm_j=vector<float>(3,0.0);
                tm_i1 =0;
                atm_idy=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
            bond0max=-1000.0;bond1max=-1000.0;bond2max=-1000.0;
            bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;

            for(int j=0;j<tmp_mat.size();j++)
            {
                if(tmp_mat[j].elt_ == "CA")
                {
                    vector<float> cartX1;
                    vector<float> fracX1; 
                //    elt_i = tmp_mat[3*tm_i1+1].elt_;
                    elt_i = tmp_mat[j].elt_;
                    elt_i = elt_i[0];
                    OneGaussianScattering sig_j = get_A( elt_i );
                    k = sig_j.k( theDensityMap.effectiveB );
                    C = sig_j.C( k );
                    if ( C < 1e-6 ) continue;  

            //        cartX1 = tmp_mat[3*tm_i1+1].x_;
                    cartX1 = tmp_mat[j].x_;
            //        tm_i1 = tm_i1 + 1;
                    MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                    atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);                 
                    for(int z=1;z<=theDensityMap.density.u3();z++)
                    {
                        atm_j[2] = z;
                        del_ij[2] =(atm_idy[tm_i1][2]-atm_j[2])/theDensityMap.grid[2];
                        if ( del_ij[2] > 0.5 ) del_ij[2]-=1.0;
                        if ( del_ij[2] < -0.5 ) del_ij[2]+=1.0;                    
                        del_ij[0] = del_ij[1] = 0.0;
                        vector<float> frac_tmpz;
                        MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpz);
                        if(square_len(frac_tmpz)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;               
                        for(int y=1;y<=theDensityMap.density.u2();y++)
                        {
                            atm_j[1] = y;
                            del_ij[1] = (atm_idy[tm_i1][1] - atm_j[1])/theDensityMap.grid[1];
                            // wrap-around??
                            if ( del_ij[1] > 0.5 ) del_ij[1]-=1.0;
                            if ( del_ij[1] < -0.5 ) del_ij[1]+=1.0;                        
                            del_ij[0] = 0.0;
                            vector<float> frac_tmpy;
                            MatrixTimesTransVector(theDensityMap.f2c,del_ij,frac_tmpy);
                            if(square_len(frac_tmpy)> (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;                      
                            for(int x=1;x<=theDensityMap.density.u1();x++)
                            {
                                atm_j[0] = x;
                                del_ij[0] = (atm_idy[tm_i1][0] - atm_j[0])/theDensityMap.grid[0];
                                // wrap-around??
                                if ( del_ij[0] > 0.5 ) del_ij[0]-=1.0;
                                if ( del_ij[0] < -0.5 ) del_ij[0]+=1.0;                            
                                vector<float> cart_del_ij2;
                                MatrixTimesTransVector(theDensityMap.f2c,del_ij,cart_del_ij2);
                                float d2 = square_len(cart_del_ij2);
                                if(d2 > (theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK_PADDING+theDensityMap.ATOM_MASK) ) continue;
                            
                                float atm = C*exp(-k*d2);
                                float sigmoid_msk = exp( d2 - (theDensityMap.ATOM_MASK)*(theDensityMap.ATOM_MASK)  );
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.ATOM_MASK));
                            //    float sigmoid_msk = exp( sqrt(d2) - (theDensityMap.CA_MASK));
                                float inv_msk = 1.0/(1.0+sigmoid_msk);
                //                rhoC2(x,y,z) += atm;
                //                rhoC2(x,y,z) += atm;
                                rhoC2(x,y,z) += atm;
                                inv_rho_mask2(x,y,z) *= (1.0 - inv_msk);


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
                    tm_i1 = tm_i1 + 1;
                }
            }
    //        rhoC2 = rhoC01;
            max_rhc = 0.0;
            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                if(rhoC2[x]>max_rhc) max_rhc = rhoC2[x] ;
            }     
    /*        for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
        //        cout<<theDensityMap.density[x]<<" ";
                if(rhoC2[x]<(3.0/5.0)*max_rhc) rhoC2[x]=0; 
            } */
//            rhoC2 = rhoC01;                           
            sumC_i2=0.0, sumO_i2=0.0, sumCO_i2=0.0, vol_i2=0.0, CC_i2=0.0;
            sumO2_i2=0.0, sumC2_i2=0.0, varC_i2=0.0, varO_i2=0.0;
            clc_x2=0.0, obs_x2=0.0, eps_x2=0.0;
    /*            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
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
                }   */
            for(int x=bond0min;x<=bond0max;x++)
            {
               for(int y=bond1min;y<=bond1max;y++)
                {
                    for(int z=bond2min;z<=bond2max;z++)
                    {
           //            tmp_den=0.0;
            //            tmp_rhc =0.0;
            //            if(theDensityMap.density[x] > (2.0/3.0)*max_den) tmp_den= theDensityMap.density[x];
            //            if(rhoC2[x]>(2.0/3.0)*max_rhc) tmp_rhc=rhoC2[x];
            //            if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        clc_x2 = rhoC2(x,y,z);
            //            clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density(x,y,z);
            //            obs_x2 = tmp_den;
                        eps_x2 = 1-inv_rho_mask2(x,y,z); //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
            //            eps_x2 = 1.0;

                        // SMOOTHED
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
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }             
    /*        sumis=0.0;
            sumu=0.0;
            for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                tmp_den = 0.0; 
                tmp_rhc=0.0;
                if(theDensityMap.density[x]>(3.0/5.0)*max_den) tmp_den=1.0;
                if(rhoC2[x]>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
//                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                // fetch this point
                sumis=sumis+tmp_den*tmp_rhc;
                sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
            }
            if ( sumu == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = sumis/sumu;
            }   */
    /*        sumis=0.0;
            sumu=0.0;
            for(int x=bond0min;x<=bond0max;x++)
            {
               for(int y=bond1min;y<=bond1max;y++)
                {
                    for(int z=bond2min;z<=bond2max;z++)
                    {
                        tmp_den = 0.0; 
                        tmp_rhc=0.0;
                        if(theDensityMap.density(x,y,z)>(3.0/5.0)*max_den) tmp_den=1.0;
                        if(rhoC2(x,y,z)>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
        //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        sumis=sumis+tmp_den*tmp_rhc;
                        sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
                    }
                }
            }          
            if ( sumu == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = sumis/sumu;
            }             */
/*            sumis=0.0;
            sumu=0.0;   
            bond0max=0.0;bond1max=0.0;bond2max=0.0;
            bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;
            tm_i1=0;
            atm_idy=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
            for(int j=0;j<tmp_mat.size();j++)
            {
 //               tmp_mat.push_back(pointsBx[j]);
                if(tmp_mat[j].elt_ == "CA")
                {
                    vector<float> cartX1;
                    vector<float> fracX1; 

             //       cartX1 = tmp_mat[3*tm_i1+1].x_;
                    cartX1 = tmp_mat[j].x_;
                    
                    MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                    atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);
                    if(atm_idy[tm_i1][0]>bond0max) bond0max = atm_idy[tm_i1][0];
                    if(atm_idy[tm_i1][1]>bond1max) bond1max = atm_idy[tm_i1][1];
                    if(atm_idy[tm_i1][2]>bond2max) bond2max = atm_idy[tm_i1][2];
                    if(atm_idy[tm_i1][0]<bond0min) bond0min = atm_idy[tm_i1][0];
                    if(atm_idy[tm_i1][1]<bond1min) bond1min = atm_idy[tm_i1][1];
                    if(atm_idy[tm_i1][2]<bond2min) bond2min = atm_idy[tm_i1][2];
                    tm_i1 = tm_i1 + 1;
                }
            }  
            for(int x=bond0min;x<=bond0max;x++)
            {
               for(int y=bond1min;y<=bond1max;y++)
                {
                    for(int z=bond2min;z<=bond2max;z++)
                    {
                        tmp_den = 0.0; 
                        tmp_rhc=0.0;
                        if(theDensityMap.density(x,y,z)>(3.0/5.0)*max_den) tmp_den=1.0;
                        if(rhoC2(x,y,z)>(3.0/5.0)*max_rhc) tmp_rhc=1.0;
        //                if(theDensityMap.density[x]<(3.0/5.0)*max_den) continue;
                        // fetch this point
                        sumis=sumis+tmp_den*tmp_rhc;
                        sumu=sumu+tmp_den+tmp_rhc-tmp_den*tmp_rhc;
                    }
                }
            }          
            if ( sumu == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = sumis/sumu;
            }      */                         
//            vector<vector<float> > bondx(3,vector<float>(2,0.0));
/*            bond0max=0.0;bond1max=0.0;bond2max=0.0;
            bond0min=1000.0;bond1min=1000.0;bond2min=1000.0;
            tm_i1=0;
            atm_idy=vector<vector<float> > ((rand_b-rand_a+1),vector<float>(3,0.0));
            for(int j=0;j<tmp_mat.size();j++)
            {
 //               tmp_mat.push_back(pointsBx[j]);
                if(tmp_mat[j].elt_ == "CA")
                {
                    vector<float> cartX1;
                    vector<float> fracX1; 

             //       cartX1 = tmp_mat[3*tm_i1+1].x_;
                    cartX1 = tmp_mat[j].x_;
                    
                    MatrixTimesTransVector(theDensityMap.c2f,cartX1,fracX1);
                    atm_idy[tm_i1][0] = pos_mod (double(fracX1[0]*theDensityMap.grid[0] - theDensityMap.origin[0] + 1) , (double)theDensityMap.grid[0]);
                    atm_idy[tm_i1][1] = pos_mod (double(fracX1[1]*theDensityMap.grid[1] - theDensityMap.origin[1] + 1) , (double)theDensityMap.grid[1]);
                    atm_idy[tm_i1][2] = pos_mod (double(fracX1[2]*theDensityMap.grid[2] - theDensityMap.origin[2] + 1) , (double)theDensityMap.grid[2]);
                    if(atm_idy[tm_i1][0]>bond0max) bond0max = atm_idy[tm_i1][0];
                    if(atm_idy[tm_i1][1]>bond1max) bond1max = atm_idy[tm_i1][1];
                    if(atm_idy[tm_i1][2]>bond2max) bond2max = atm_idy[tm_i1][2];
                    if(atm_idy[tm_i1][0]<bond0min) bond0min = atm_idy[tm_i1][0];
                    if(atm_idy[tm_i1][1]<bond1min) bond1min = atm_idy[tm_i1][1];
                    if(atm_idy[tm_i1][2]<bond2min) bond2min = atm_idy[tm_i1][2];
                    tm_i1 = tm_i1 + 1;
                }
            } 
            for(int x=bond0min;x<=bond0max;x++)
            {
               for(int y=bond1min;y<=bond1max;y++)
                {
                    for(int z=bond2min;z<=bond2max;z++)
                    {
                        clc_x2 = rhoC2(x,y,z);
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density(x,y,z);
        //                obs_x2 = tmp_den;
        //                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                        eps_x2 = 1.0;

                        // SMOOTHED
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
            varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2 / vol_i2 ) ;
            if ( varC_i2 == 0 || varO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2/ vol_i2) / sqrt( varC_i2 * varO_i2 );
            }                */
 /*               for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x ) {
                    clc_x2 = rhoC2[x];
                    obs_x2 = theDensityMap.density[x];
                    eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal

                    // SMOOTHED
                    sumCO_i2 += eps_x2*clc_x2*obs_x2;
                    sumO_i2  += eps_x2*obs_x2;
                    sumO2_i2 += eps_x2*obs_x2*obs_x2;
                    sumC_i2  += eps_x2*clc_x2;
                    sumC2_i2 += eps_x2*clc_x2*clc_x2;
                    vol_i2   += eps_x2;
                }
                sumO_i2 = sumO_i2/vol_i2;
                sumC_i2 = sumC_i2/vol_i2;
                sumCO_i2 = sumCO_i2/vol_i2;
                sumC2_i2 = sumC2_i2/vol_i2;
                sumO2_i2 = sumO2_i2/vol_i2;
                varC_i2 = (sumC2_i2 - sumC_i2*sumC_i2);
                varO_i2 = (sumO2_i2 - sumO_i2*sumO_i2) ;
                if ( varC_i2 == 0 || varO_i2 == 0 ) {
                    CC_i2 = 0;
                } else {
                    CC_i2 = (sumCO_i2 - sumC_i2*sumO_i2) / sqrt( varC_i2 * varO_i2 );
                }  */
    /*        for(int x=bond0min;x<=bond0max;x++)
            {
               for(int y=bond1min;y<=bond1max;y++)
                {
                    for(int z=bond2min;z<=bond2max;z++)
                    {
                        clc_x2 = rhoC2(x,y,z);
        //                clc_x2 = tmp_rhc;
                        obs_x2 = theDensityMap.density(x,y,z);
        //                obs_x2 = tmp_den;
        //                eps_x2 = 1-inv_rho_mask2[x]; //1/(1+exp( (0.01-rho_calc(x,y,z)) * 1000 ));  // sigmoidal
                        eps_x2 = 1.0;

                        // SMOOTHED
                        vol_i2 += eps_x2;
                        sumC_i2 += (eps_x2*clc_x2*clc_x2);
                        sumO_i2 += (eps_x2*obs_x2*obs_x2);
                        sumCO_i2 += (eps_x2*clc_x2*obs_x2);
                    }
                }
            }          
            if ( sumC_i2 == 0 || sumO_i2 == 0 ) {
                CC_i2 = 0;
            } else {
                CC_i2 = sumCO_i2/(sqrt(sumC_i2+sumO_i2)); 
            }             */
/*                for ( int x=0; x<theDensityMap.density.u1()*theDensityMap.density.u2()*theDensityMap.density.u3(); ++x )
                {
                    eps_x2 = 1-inv_rho_mask2[x];
                    vol_i2 += eps_x2;
                    sumC_i2 += (eps_x2*rhoC2[x]*rhoC2[x]);
                    sumO_i2 += (eps_x2*theDensityMap.density[x]*theDensityMap.density[x]);
                    sumCO_i2 += (eps_x2*rhoC2[x]*theDensityMap.density[x]);
                }
            //    sumCO_i2 = sumCO_i2/vol_i2;
            //    sumO_i2 = sumO_i2/vol_i2;
            //    sumC_i2 = sumC_i2/vol_i2;
                CC_i2 = sumCO_i2/(sqrt(sumC_i2*sumO_i2));  */

//                cout<< "DDDDDXX"<<endl;
                // ss inforamtion     
                float new_Ehbond=0.0;
        //        point3f *new_decstr;
                numseq=abs(rand_b-rand_a+1);
                vector<point3f> new_decstr(numseq);

            /*    vector<boneinfo> bby(numseq);
                tmp_tt=0;
                for(int j=rand_a;j<=rand_b;j++)
                {
                    bby[tmp_tt]=bb[j];
                    tmp_tt = tmp_tt +1;
                }       */  
                if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
                {                       
                /*    tmp_tt=0;
                    for(int j=0;j<numseq;j++)
                    {
                        bby[tmp_tt].indn=3*tmp_tt+0;
                        bby[tmp_tt].indca=3*tmp_tt+1;
                        bby[tmp_tt].indc=3*tmp_tt+2;   
                        tmp_tt=tmp_tt+1;
                    }     
                    calcsse2(bby, numseq, tmp_mat);    */

//                if(numtoss(bb[randtx].sst)=='H'||numtoss(bb[randtx].sst)=='E')
//                {
                    tmp_tt=0;
                    for(int j=rand_a;j<=rand_b;j++)
                    {
                    //    new_decstr[tmp_tt].ss2=numtoss(bby[j].sst);
                        new_decstr[tmp_tt].ss2=numtoss(bb[j].sst);
                        new_decstr[tmp_tt].x=tmp_mat[3*tmp_tt+1].x_[0];
                        new_decstr[tmp_tt].y=tmp_mat[3*tmp_tt+1].x_[1];
                        new_decstr[tmp_tt].z=tmp_mat[3*tmp_tt+1].x_[2];
                        new_decstr[tmp_tt].ptn.x=tmp_mat[3*tmp_tt+0].x_[0];
                        new_decstr[tmp_tt].ptn.y=tmp_mat[3*tmp_tt+0].x_[1];
                        new_decstr[tmp_tt].ptn.z=tmp_mat[3*tmp_tt+0].x_[2];
                        new_decstr[tmp_tt].ptc.x=tmp_mat[3*tmp_tt+2].x_[0];
                        new_decstr[tmp_tt].ptc.y=tmp_mat[3*tmp_tt+2].x_[1];
                        new_decstr[tmp_tt].ptc.z=tmp_mat[3*tmp_tt+2].x_[2];
                        tmp_tt=tmp_tt+1;
                    } 
                    new_Ehbond=energyhbondcanc(new_decstr,numseq);
                }            
            //    cout<<"new_Ehbond: "<<new_Ehbond<<endl; 
                tm_i1=0;                    
    /*            float new_ss_E =0.0;
                vector<string > new_ss_inf(CA_num);
                vector<vector<float> > new_CA_cor;
                tm_i1=0;                
                for(int j=0;j<tmp_mat.size();j++)
                {
                    if(tmp_mat[j].elt_ == "CA")
                    {
                        new_CA_cor.push_back(tmp_mat[j].x_);
                        tm_i1 = tm_i1+1 ; 
                    }                 
                }
                make_secx(new_CA_cor,CA_num,new_ss_inf);
                ss_ab_num=0; // helix and b strand number
                ss_c_num =0; // coil number                  
                for(int j=0;j<new_ss_inf.size();j++)
                {
                    if(new_ss_inf[j] == 'H')
                    {
                        new_ss_E = new_ss_E + 0.1;
                        ss_ab_num = ss_ab_num +1;
                    }
                    else if(new_ss_inf[j] == 'E')
                    {
                        new_ss_E = new_ss_E + 0.1;
                        ss_ab_num = ss_ab_num +1;
                    }
                    else{
                        new_ss_E = new_ss_E + 0.2;
                        ss_c_num = ss_c_num +1;
                    }                    
                }
                new_ss_E = new_ss_E/(CA_num);    */             
                new_CC =1.0-CC_i2; 
        //        cout<< "new_CC: "<<new_CC<<" ";
            //    new_dE = 30*new_CC + new_vwd/((pnumx-1)*(rand_b-rand_a+1));
            //    new_dE = 0.8*new_CC + 0.5*new_clash/abs(rand_b-rand_a+1) + 0.5*new_Ehbond/abs(rand_b-rand_a+1);
                new_dE = 0.5*new_CC + 0.5*new_clash/abs(rand_b-rand_a+1) + 0.5*new_Ehbond/abs(rand_b-rand_a+1);
            //    cout<<"new coor:  "<<tmp_mat[5].x_[0]<<" "<<tmp_mat[5].x_[1]<<" "<<tmp_mat[5].x_[2]<<endl;
                cout<<"old: "<< old_CC<<" " <<old_clash<<" "<<old_Ehbond<<endl;
                cout<<"new: "<< new_CC<<" " <<new_clash<<" "<<new_Ehbond<<endl;
                dE = new_dE - old_dE;
//                dE = 1.0/new_dE - 1.0/old_dE;
                if(new_dE <old_dE)
                {
                    int tmp_tmx=0;
                    for(int j=rand_a;j<=rand_b;j++)
                    {
                        pointsBx[3*j+0] = tmp_mat[3*tmp_tmx+0];
                        pointsBx[3*j+1] = tmp_mat[3*tmp_tmx+1];
                        pointsBx[3*j+2] = tmp_mat[3*tmp_tmx+2];
                    //    bb[j]=bby[tmp_tmx];
                        tmp_tmx = tmp_tmx +1;
                    }
                //    rhoC0 =rhoC01;
                    cout<<" new_dE: "<<new_dE<<endl; 
                }                
                else
                {
                //    float tmpx=rand()/double(RAND_MAX);
                    float tmpx=randf0and1();
                    float mc_v = exp(-dE/(KT));
                    if(tmpx < mc_v)
                    {
                        int tmp_tmx=0;
                        for(int j=rand_a;j<=rand_b;j++)
                        {
                            pointsBx[3*j+0] = tmp_mat[3*tmp_tmx+0];
                            pointsBx[3*j+1] = tmp_mat[3*tmp_tmx+1];
                            pointsBx[3*j+2] = tmp_mat[3*tmp_tmx+2];
                        //    bb[j]=bby[tmp_tmx];
                            tmp_tmx = tmp_tmx +1;
                        }
                    //    rhoC0 =rhoC01;
                        cout<<" XXX new_dE: "<<new_dE<<endl; 
                    }       
                }  
//                tmp_mat.clear(); 

    //    }

    } // the recile 
}

/*    vector<int > beg_end;
    getLinLen(randop,pointsBx,bb,beg_end);
    cout<<"beg end: "<<beg_end[0]<<" "<<beg_end[1]<<endl;

    vector<poseCoord> fra_tmp;
    for(int tt=beg_end[0];tt<=beg_end[1];tt++)
    {
        fra_tmp.push_back(pointsBx[3*tt+0]);
        fra_tmp.push_back(pointsBx[3*tt+1]);
        fra_tmp.push_back(pointsBx[3*tt+2]);
    }
    connectgap(fra_tmp);
    int tmpee=0;
    for(int tt=beg_end[0];tt<=beg_end[1];tt++)
    {
        pointsBx[3*tt+0]=fra_tmp[3*tmpee+0];
        pointsBx[3*tt+1]=fra_tmp[3*tmpee+1];
        pointsBx[3*tt+2]=fra_tmp[3*tmpee+2];
        tmpee = tmpee +1;
    }   */ 

    
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
    string inputMRC;
//  poseCoords pose1,pose2;
    Model pose;

    clock_t start,finish;
    double totaltime;
    start=clock();  

//  readPDB
    strcpy(inputPDB1,argv[1]);
    readpdbstructurex(inputPDB1,pose);

    poseCoords posey;
//    readPDBcoords(inputPDB1,posey);
    cout<< "Read PDB file : "<<inputPDB1<<endl;

    for(int i=0;i<5;i++)
    {
        for(int j=0;j<5;j++)
        {
            cout<<randIntCustom(0,20)<<" "<<randf0and1()<<endl;
        }
    }


    float MRC_reso;
    float mapsampling = 0.0;
    inputMRC = argv[2];
    MRC_reso = atof(argv[3]);
    mapsampling = atof(argv[4]);
    int dm_beg= atoi(argv[5]);
    int dm_end = atoi(argv[6]);
    char bindir[600];
    strcpy(bindir,argv[7]);

//  strcpy(inputPDB2,argv[2]);
//  readPDBcoordss(inputPDB1,pose1);
//  cout<<"000  "<<endl;
//  readPDBcoordss(inputPDB2,pose2);
//  readpdbstructurex(inputPDB2,pose2);
//  cout<< "11 "<<endl;
//  cout<<" 11 "<<endl;
    Model fin_p,fin_px,fin_py;  
    cout<<" 11 "<<endl;
    fin_p=REMC_sampley(pose,inputMRC,MRC_reso,mapsampling,dm_beg,dm_end);
/*    string outnamex="11111.pdb";
    writePDBStructure(inputPDB1,fin_p,outnamex);
    cout<<"0000"<<endl;
    char cmd[1000];
    sprintf(cmd, "./Scwrl4 -i 11111.pdb -o 22222.pdb");
    system(cmd);

    cout<<"2222"<<endl;

    char rmcmd[100];
    sprintf(rmcmd, "rm -f 11111.pdb");
    system(rmcmd);

    cout<<"3333"<<endl;

    string inputPDB2="22222.pdb";
    Model poseyy;
    readpdbstructurex(inputPDB2.data(),poseyy);   
    sprintf(rmcmd, "rm -f 22222.pdb");
    system(rmcmd); */

//    string outname;
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

//    string outnamex="11111.pdb";
    writePDBStructure(inputPDB1,fin_p,outnamex);
    cout<<"0000"<<endl;
    char cmd[1000];
    sprintf(cmd, "%s/pulchra -e %s",bindir,outnamex.c_str());
    system(cmd);

    cout<<"2222"<<endl;

    char rmcmd[100];
    sprintf(rmcmd, "rm -f %s",outnamex.c_str());
    system(rmcmd);

    cout<<"3333"<<endl;

    string inputPDB2=tmp_rx+ "_11.rebuilt.pdb";
    Model poseyy;
    readpdbstructurex(inputPDB2.data(),poseyy);   
    sprintf(rmcmd, "rm -f %s",inputPDB2.c_str());
    system(rmcmd); 
    


    cout<< "444"<<endl;
    fin_py = poseyy;  
//  fin_px=REMC_samplex(fin_p,pose2);
//  fin_py=REMC_samplexy(fin_px,pose2);
//  writePDBcoordss(inputPDB1, fin_p);
    string outname;
/*    vector<string> tmp_str,tmp_strx;
    tmp_str = string_splitx(argv[1],'.');
    tmp_strx = string_splitx(tmp_str[0],'/'); */

    if(tmp_strx.size()>0)
    {
        outname = tmp_strx[tmp_strx.size()-1] + "_outp00x03.pdb";
    }
    else
    {
        outname = tmp_strx[0] + "_outp00x03.pdb";
    }
    
    writePDBStructure(inputPDB1,fin_py,outname);
    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"\nrunning time: "<<totaltime<<" s！"<<endl;  
    return 0;
}