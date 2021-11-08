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
//      cout<<nn << " ";
        if(buf.length()<54) continue;
//      if ( buf.substr(0,4) =="ATOM" && buf.substr(12,4) ==" CA " ) {
        if ( buf.substr(0,4) =="ATOM" && buf.substr(13,1) !="H" && buf.substr(13,3) !="OXT") {
            poseCoord atom_i;

            atom_i.x_.push_back(atof(buf.substr(30,8).c_str()));
            atom_i.x_.push_back(atof(buf.substr(38,8).c_str()));
            atom_i.x_.push_back(atof(buf.substr(46,8).c_str()));
            atom_i.B_ = 0.0;// atof(buf.substr(60,6).c_str())
            atom_i.elt_ = buf.substr(12,4);

//          nn=nn+1;
//          cout<< nn<<" ";
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

/*  for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<A[i].size();j++)
        {
//          dis = dis + (A[i][j]-B[i][j])*(A[i][j]-B[i][j]);    
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
//  float F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang);

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
//  cout<<"temp: "<<temp<<endl;

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


int main(int argc, char* argv[])
{
    char inputPDB1[600]; 
//  poseCoords pose1,pose2;
    Model posex;

//  readPDB
    strcpy(inputPDB1,argv[1]);
    readpdbstructurex(inputPDB1,posex);

    cout<< "Read PDB file : "<<inputPDB1<<endl;

    string outnamex;
    vector<string> tmp_str,tmp_strx;
    tmp_str = string_splitx(argv[1],'.');
    tmp_strx = string_splitx(tmp_str[0],'/'); 
    string tmp_rx;
    if(tmp_strx.size()>0)
    {
        outnamex = tmp_strx[tmp_strx.size()-1];
        tmp_rx = tmp_strx[tmp_strx.size()-1];
    }
    else
    {
        outnamex = tmp_strx[0] ;
        tmp_rx = tmp_strx[0];
    }


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

            }

            pointsx.push_back(rotx);
//          pointsy.push_back(roty);
            pnum = pnum +1; 
//          cout<< pnum<<endl;      
        }
    }
    cout<<"pnum,pnumx: "<<pnum<<" "<<pnumx<<endl;
    pointsBx0 = pointsBx;  

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

    int ppxxx=0;
    if(ppxxx == 0)
    {
        int numseq=pnum;
        vector<point3f> decstr(numseq);                
        int tmp_tt=0;
        for(int j=0;j<pnum;j++)
        {
            decstr[tmp_tt].ss2 = numtoss(bb[j].sst);// predicted secondary structure
            decstr[tmp_tt].ssm = numtoss(bb[j].sst);
            decstr[tmp_tt].x = pointsBx[3*tmp_tt+1].x_[0];
            decstr[tmp_tt].y = pointsBx[3*tmp_tt+1].x_[1];
            decstr[tmp_tt].z = pointsBx[3*tmp_tt+1].x_[2];
            decstr[tmp_tt].ptn.x = pointsBx[3*tmp_tt+0].x_[0];
            decstr[tmp_tt].ptn.y = pointsBx[3*tmp_tt+0].x_[1];
            decstr[tmp_tt].ptn.z = pointsBx[3*tmp_tt+0].x_[2];
            decstr[tmp_tt].ptc.x = pointsBx[3*tmp_tt+2].x_[0];
            decstr[tmp_tt].ptc.y = pointsBx[3*tmp_tt+2].x_[1];
            decstr[tmp_tt].ptc.z = pointsBx[3*tmp_tt+2].x_[2];
            tmp_tt=tmp_tt+1;
        } 
    //    old_Ehbond=energyhbondcanc(decstr,numseq);
        vector<sssegment> sse;
        int numsse=0;
        genesse(decstr,numseq,sse,numsse) ;
        ofstream outfile;
        string fin_name = outnamex + ".sse"; 
        outfile.open(fin_name.c_str());
        if(!outfile.is_open())
        {
            cout<<"error to open the sse file"<<endl;
        }
        for(int j=0;j<numsse;j++)
        {
            outfile<<sse[j].ss<<"   "<<sse[j].init<<"      "<<sse[j].term<<endl;
        }
        outfile.close();
    } 

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

    return 0;
}