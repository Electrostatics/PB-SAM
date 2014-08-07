#include <iostream>
#include <fstream>
#include <vector>
#include "util_sphere.h"
#include "string.h"

// various functions to read files


// Reads in the coordinates of a set of points that describe the surface of 
// the molecule. Also reads in the normals to the surface at each point.
// The data is read from a ".vert" file output by the program MSMS
void 
readSurface(const char * fname, vector<CPnt> & SP, vector<CPnt> & NP )
{
  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open surface file " << fname << endl;
      exit(0);
    }

  CPnt scen;

  string s;
  getline(fin, s);
  while (s[0] == '#')
    getline(fin, s);

  REAL x,y,z,nx,ny,nz;
  getline(fin, s);
  while (!fin.eof())
    {
      istringstream ss(s);
      
      ss >> x >> y >> z;
      SP.push_back(CPnt(x,y,z));
      scen += CPnt(x,y,z);

      ss >> nx >> ny >> nz;
      NP.push_back(CPnt(nx,ny,nz));
      getline(fin, s);
    }

  scen /= SP.size();
  cout<< SP.size() << " surface points read!!" << endl;
  cout<< "estimated center from mol surface: "<<scen<<endl;
}

// Read an PQR format file
void readpqr(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
	     vector<REAL> & rad, vector<CPnt> & cen)
{
  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open file " << fname << endl;
      exit(0);
    }

  pnt.resize(0);
  ch.resize(0);
  rad.resize(0);
  cen.resize(0);

  char buf[100], temp[10];
  REAL sum = 0.0;
  vector<REAL> R;
  pnt.clear();
  ch.clear();
  rad.clear();
  cen.clear();
  fin.getline(buf,99);

  int iCen, iCoord, iCharge, iRad;
  iCen = 17;
  iCoord = 31;
  iCharge = 55;
  iRad = 63;
  /*
  iCen = 18;
  iCoord = 31;
  iCharge = 57;
  iRad = 65;
  */ 

  while (!fin.eof())
    {
      double x,y,z,c,r;
      if (strncmp(&(buf[0]),"ATOM",4) == 0)
	{
	  sscanf(&(buf[iCoord]), "%lf %lf %lf", &x, &y, &z); // position
	  strncpy(temp, &(buf[iRad]), 6); // radius
	  temp[6] = 0;
	  if (strcmp(temp, "      ") == 0)
	    r = 0.0;
	  else
	    sscanf(temp, "%lf", &r);

	  sscanf(&(buf[iCharge]), "%lf", &c); // charge
	  
	  // read in as centers that specifies dielectric boundary
	  if (strncmp(&(buf[iCen]),"CEN",3) == 0)
	    {
	      rad.push_back(r);
	      cen.push_back(CPnt(x,y,z));
	    }

	  // read in as atoms
	  else
	    {
	      pnt.push_back(CPnt(x,y,z));
	      ch.push_back(c);
	      R.push_back(r);
	      sum += c;
	    }
	}

      fin.getline(buf,99);     
    }

  // for(int i=0; i<pnt.size(); i++) cout <<"charge "<<i+1<<" "<<pnt[i]<<" "<<ch[i]<<endl;
  for(int i=0; i<cen.size(); i++) cout <<"cen "<<i+1<<" "<<cen[i]<<" "<<rad[i]<<endl;

  cout << fname << ": atoms: " << ch.size() << ", net charge: " 
       << sum <<endl;
  for(int k = 0; k < rad.size(); k++) cout <<k<<" boundary sradius: " << rad[k] << ", scen: " << cen[k] << endl; 
}

// Read an PQR format file
void readpqrAtoms(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
	     vector<REAL> & R)
{
  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open file " << fname << endl;
      exit(0);
    }

  pnt.resize(0);
  ch.resize(0);
  R.resize(0);

  char buf[100], temp[10];
  int nc = 0;
  REAL sum = 0.0;

  pnt.clear();
  ch.clear();
  R.clear();

  fin.getline(buf,99);

  int iCen, iCoord, iCharge, iRad;
   
  iCen = 18;
  iCoord = 31;
  iCharge = 57;
  iRad = 65;

  /*iCen = 17;
  iCoord = 31;
  iCharge = 55;
  iRad = 63;
	*/ 
  while (!fin.eof())
    {
      double x,y,z,c,r;
      if (strncmp(&(buf[0]),"ATOM",4) == 0)
	{
	  sscanf(&(buf[iCoord]), "%lf %lf %lf", &x, &y, &z); // position
	  sscanf(&(buf[iCharge]), "%lf", &c); // charge
	  sscanf(&(buf[iRad]), "%lf", &r);
	  
	  pnt.push_back(CPnt(x,y,z));
	  ch.push_back(c);
	  R.push_back(r);
	  sum += c;
	  if(fabs(c) > 1e-15) nc++;
 	}

      fin.getline(buf,99);     
    }

  // for(int i=0; i<pnt.size(); i++) cout <<"charge "<<i+1<<" "<<pnt[i]<<" "<<ch[i]<<endl;
  cout << fname << ": Total atoms: " << pnt.size()<<" charged atoms: "<<nc<< ", net charge: " 
       << sum <<endl;
}

// Read an MDP format file
void readmdp(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
	     vector<REAL> & rad, vector<CPnt> & cen)
{
  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open file " << fname << endl;
      exit(0);
    }

  char buf[100], temp[10];
  REAL sum = 0.0;
  vector<REAL> R;
  pnt.clear();
  ch.clear();
  rad.clear();
  cen.clear();
  fin.getline(buf,99);

  while (!fin.eof())
    {
      double x,y,z,c,r;
      if (strncmp(&(buf[0]),"ATOM",4) == 0)
	{
	  sscanf(&(buf[31]), "%lf %lf %lf", &x, &y, &z);
	  strncpy(temp, &(buf[55]), 6);
	  temp[6] = 0;
	  if (strcmp(temp, "      ") == 0)
	    r = 0.0;
	  else
	    sscanf(temp, "%lf", &r);

	  sscanf(&(buf[61]), "%lf", &c);
	  
	  // specifies dielectric boundary
	  if (strncmp(&(buf[17]),"DUM",3) == 0)
	    {
	      rad.push_back(r);
	      cen.push_back(CPnt(x,y,z));
	    }
	  else
	    {
	      pnt.push_back(CPnt(x,y,z));
	      ch.push_back(c);
	      R.push_back(r);
	      sum += c;
	    }
	}

      fin.getline(buf,99);     
    }

  REAL w = 0.0;
  CPnt S;
  for (int i = 0; i < pnt.size(); i++)
    {
      REAL wt = R[i]*R[i]*R[i];
      S += wt*pnt[i];
      w += wt;
    }

  S *= (1/w);

  REAL mx = 0.0;
  for (int i = 0; i < pnt.size(); i++)
    {
      REAL d = (pnt[i] - S).norm();
      if (ch[i] != 0.0 && d > mx)
	mx = d;
    }

  // cout << S << " " << mx << endl;
  cout << fname << ": atoms: " << ch.size() << ", net charge: " 
       << sum <<endl;
  for(int k = 0; k < rad.size(); k++) cout <<k<<" boundary sradius: " << rad[k] << ", scen: " << cen[k] << endl; 
}
// Write out an MDP file 
// It is just like a pdb file only the last two entries after the 
// X,Y,Z coordinates are the radius and the charge of the atom 
void writemdp(const char * ofname, const vector<CPnt> & scen, double * srad,
	      const vector<CPnt> & pnt, const vector<double> ch)
{
  ofstream fout(ofname);

  if (!fout.is_open())
    {
      cout << "Could not open file " << ofname << endl;
      exit(0);
    }

  char buf[50];
  for (int i = 0; i < pnt.size(); i++)
    {    
      fout << "ATOM  ";

      fout.width(5);
      fout.setf(ios::right, ios::adjustfield);
      fout << i << ' ';

      fout.setf(ios::left, ios::adjustfield);
      char aname[5];
      if (ch[i] > 0)
	 strcpy(aname, " C1 ");
      else
	strcpy(aname, " C2 ");
      fout << aname << ' ';
      fout << "CAV" << ' ' << "A";

      fout.setf(ios::right, ios::adjustfield);
      fout.width(4);
      fout << 1 << "    ";
  
      sprintf(buf, "%8.3f%8.3f%8.3f%8.3f%8.3f", pnt[i].x(), pnt[i].y(), 
	      pnt[i].z(), 0.0, ch[i]);
      fout << buf << endl;;
    }

  for (int i = 0; i < scen.size(); i++)
    {    
      fout << "ATOM  ";

      fout.width(5);
      fout.setf(ios::right, ios::adjustfield);
      fout << i+ch.size() << ' ';

      fout.setf(ios::left, ios::adjustfield);
      char aname[5] = " O  ";
      fout << aname << ' ';
      fout << "CAV" << ' ' << "A";

      fout.setf(ios::right, ios::adjustfield);
      fout.width(4);
      fout << 1 << "    ";
  
      sprintf(buf, "%8.3f%8.3f%8.3f%8.3f%8.3f", scen[i].x(), scen[i].y(), 
	      scen[i].z(), srad[i],0.0);
      fout << buf << endl;;
    }
}

// Print out the entries of a matrix
void printM(double * A, int m, int n)
{
  for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
	cout << A[j*m+i] << " ";
      cout << endl;
    }
  cout << endl;
}
