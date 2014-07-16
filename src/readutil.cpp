#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "string.h"
#include "readutil.h"
#include "util.h"

//!  readgrid function 
/*! A function to Read a DelPhi grid file (The length of 
each side of the 3D grid is GSIZE)
	\param fname a character string to read in grid from. 
	\param V a pointer to a vector of vector points
	\param scale a pointer to a vector of scaling factors
	\param pGridRef a vector of reference points */
void 
readgrid(const char * fname, float * V, float * scale, CPnt & pGridRef)
{
  FILE * fid = fopen(fname, "r");
  if (!fid)
	{
		cout << "Could not open file " << fname << endl;
		exit(0);
	}
	
  char buf[100];
  float midt[3];
	
  fread(buf,4, sizeof(char), fid);
  fread(buf, 20, sizeof(char), fid);
  fread(buf, 7, sizeof(char), fid);
  fread(buf, 10, sizeof(char), fid);
	fread(buf, 60, sizeof(char), fid);
  fread(buf, 5, sizeof(char), fid);
  fread(V, GSIZE_CB, sizeof(float), fid);
  fread(buf, 12, sizeof(char), fid);
  fread(buf, 16, sizeof(char), fid);
  fread(buf, 8, sizeof(char), fid);
  fread(scale, 1, sizeof(float), fid);
  fread(midt, 3, sizeof(float), fid);
  
  pGridRef = CPnt(midt[0],midt[1],midt[2]);
  cout << "readgrid: scale: " <<  *scale << "  pGridRef (=center for Delphi): " << pGridRef<< endl;
	
  fclose(fid);
  return;
}

//!  readSurface function 
/*! A function to read in a MSMS file. Reads in the 
		coordinates of a set of points that describe the surface of 
		the molecule. Also reads in the normals to the surface at each point.
		The data is read from a ".vert" file output by the program MSMS
	\param fname a character string to read in MSMS from. 
	\param SP a vector of xyz coordinate objects to store coordinates in
	\param NP a vector of CPnt objects to store MSMS information from */
void 
readSurface( const char * fname, vector<CPnt> & SP, vector<CPnt> & NP  )
{	
	// File handling
	ifstream fin( fname );

	// if the program is unable to open the file,
	// die
	if ( !fin.is_open( ))
	{
		cout << "Could not open surface file " << fname << endl;
		exit( 0 );
	}

	CPnt scen;
	string s;

	// Ignoring any comments in the file
	getline( fin, s );
	while ( s[0] == '#' )
		getline( fin, s );

	REAL x,y,z,nx,ny,nz;
	getline( fin, s );

	// While there is no EOF read
	while ( !fin.eof( ))
	{
		istringstream ss( s );

		ss >> x >> y >> z;
		SP.push_back( CPnt(x,y,z ));
		scen += CPnt( x,y,z );

		ss >> nx >> ny >> nz;
		NP.push_back( CPnt(nx,ny,nz ));
		getline( fin, s );
	}

	scen /= SP.size(  );
	cout<< SP.size(  ) << " surface points read!!" << endl;
	cout<< "estimated center from mol surface: "<<scen<<endl;
}	// end readSurface

//!  readpqrAtoms function 
/*! A function to read in atoms from a PQR file
	\param fname a character string to read in centers from
	\param pnt a vector of xyz coordinate objects to store coordinates in
	\param ch a vector of charges
	\param R a vector of radii of atoms */
void 
readpqrAtoms( const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
									vector<REAL> & R )
{
	
	// File handling
	ifstream fin( fname );
	
	// Die if the file wouldn't open
	if ( !fin.is_open( ))
	{
		cout << "Could not open file " << fname << endl;
		exit( 0 );
	}
	
	// Declaring variables
	pnt.resize( 0 );
	ch.resize( 0 );
	R.resize( 0 );
	
	char buf[100], temp[10];
	int nc = 0;
	REAL sum = 0.0;
	
	pnt.clear(  );
	ch.clear(  );
	R.clear(  );
	
	fin.getline( buf,99 );
	
	// Locations for coord, charge, rad
	// in the PQR file
	int iCen, iCoord, iCharge, iRad;
	
	//	iCen = 17;
	//	iCoord = 31;
	//	iCharge = 55;
	//	iRad = 63;
	
	// Other format version
	
	iCen = 18;
	iCoord = 31;
	iCharge = 58;
	iRad = 67;
	
	
	// While not at EOF
	while ( !fin.eof( ))
	{
		double x,y,z,c,r;
		// If it is an ATOM
		if ( strncmp(&(buf[0] ),"ATOM",4) == 0)
		{
			// Find the coords, charge and radius from PQR
			sscanf( &(buf[iCoord] ), "%lf %lf %lf", &x, &y, &z); // position
			sscanf( &(buf[iCharge] ), "%lf", &c); // charge
			sscanf( &(buf[iRad] ), "%lf", &r);
			
			// Add them to current list
			pnt.push_back( CPnt(x,y,z ));
			ch.push_back( c );
			R.push_back( r );
			sum += c;
			
			// Increment the number of charged atoms
			// if necessary
			if( fabs(c ) > 1e-15) nc++;
		}
		
		// Get another line
		fin.getline( buf,99 );     
	}
	
	// Printing out total charge and the total atom count etc
	cout << fname << ": Total atoms: " << pnt.size(  )
	<<" charged atoms: "<<nc<< ", net charge: " 
	<< sum <<endl;
}	// end readpqrAtoms

//!  readpqrCenters function 
/*! A function to read in CG spheres from a PQR file
	\param fname a character string to read in centers from
	\param scen a vector of xyz coordinate objects to store coordinates in
	\param srad a vector of floating points of center radii
	\param bNew a boolean of whether or not to use new formatting or not for PQR */
void readpqrCenters(const char * fname, vector<CPnt> & scen, vector<REAL> & srad, bool bNew)
{
  ifstream fin(fname);
  if (!fin.is_open())
	{
		cout << "Could not open file " << fname << endl;
		exit(0);
	}
	
  scen.clear();
  srad.clear(); 
	
  int iCoord, iCharge, iRad, iCen;
  if(bNew)
	{
		iCen = 17;
		iCoord = 27;
		iCharge = 52;
		iRad = 59;
	}
  else
	{
		iCen = 18;
		iCoord = 31;
		iCharge = 57;
		iRad = 65;
	}
	
  char buf[100], temp[10];
	
  fin.getline(buf,99);
	
  while (!fin.eof())
	{
		double x,y,z,c,r;
		if (strncmp(&(buf[0]),"ATOM",4) == 0)
		{
			
			sscanf(&(buf[iCoord]), "%lf %lf %lf", &x, &y, &z); // position
			sscanf(&(buf[iCharge]), "%lf", &c); // charge
			sscanf(&(buf[iRad]), "%lf", &r);
			
			if (strncmp(&(buf[iCen]),"CEN",3) == 0)
	    {
	      srad.push_back(r);
	      scen.push_back(CPnt(x,y,z));
	    }
			
		}
		
		fin.getline(buf,99);     
	}
	
} // end readpqrCenters

//!  readpqr function 
/*! A function to read in a PQR file that contains both atoms
and CG spheres
	\param fname a character string to read in centers from
	\param pnt a vector of xyz coordinate objects to store coordinates in
	\param ch a vector of atom partial charges
	\param rad a vector of floating points of center radii
	\param cen a vector of CG sphere xyz coordinates */
void readpqr(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
						 vector<REAL> & rad, vector<CPnt> & cen)
{
  cout <<"Reading PQR file "<<fname<<endl;
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
  int iCoord, iCharge, iRad, iCen;
	
  // center files
  iCen = 18;
  iCoord = 31;
  iCharge = 57;
  iRad = 65;
  /*
	 iCen = 17;
	 iCoord = 31;
	 iCharge = 55;
	 iRad = 62;
	 */
  while (!fin.eof())
	{
		double x,y,z,c,r;
		if (strncmp(&(buf[0]),"ATOM",4) == 0)
		{
			sscanf(&(buf[iCoord]), "%lf %lf %lf", &x, &y, &z); // position
			strncpy(temp, &(buf[iRad]), 7); // radius
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
	
  cout << fname << ": atoms: " << ch.size() << ", net charge: " 
	<< sum <<endl;
	
  REAL sphereTol = 2.0;
  CPnt rcen = computeAtomCenter(pnt);
  REAL maxR = computeMaxRadius_ByAtom(pnt, rcen, sphereTol);
	
} // end readpqr


//!  printM function 
/*! A function to print out entries of a matrix
	\param A matrix to print
	\param m number of rows of the matrix
	\param n number of columns of the matrix */
void printM( double * A, int m, int n )
{
	for ( int i = 0; i < m; i++ )
	{
		for ( int j = 0; j < n; j++ )
		{
			cout << A[j*m+i] << " ";
		}
		cout << endl;
	}
	cout << endl;
} // end printM


//!  writeMat_Unformatted function 
/*! A function to write out an unformatted matrix
	\param mat a matrix to print
	\param p number of poles of the matrix
	\param fname a file to write the matrix out to */
void writeMat_Unformatted(double *mat, int p, char* fname)
{
  ofstream fout; 
  fout.open(fname, ofstream::binary); //for writing
  if (!fout) 
	{       
		cout << "file "<< fname << " could not be opened."<< endl;
		exit(1);
	} 
	
  int length = p*p*p*p;
  fout.write( reinterpret_cast<char const *> (&p), sizeof(p)); // pole order
  fout.write( reinterpret_cast<char const *>(mat), length*sizeof(double) ); //mat
	
  fout.close();  
}

//!  readMat_Unformatted function 
/*! A function to read in a formatted matrix
	\param mat a matrix to read into
	\param maxp number of poles of the matrix
	\param fname a file to read the matrix from */
void readMat_Unformatted(double *mat, int maxp, char* fname)
{
  //  cout << "Reading unformatted imat from "<<fname<<endl;
  const int psq = maxp * maxp;
  const int pquad = psq * psq;
	
  ifstream fin;
  fin.open(fname, ifstream::binary); //for reading
  if(!fin)
	{
		cout <<"Error reading Mat:  cannot open "<<fname<<endl;
		exit(1);
	}
	
  // read pole order
  int p;
  fin.read( reinterpret_cast<char*> (&p), sizeof(int)); 
  if( p < maxp ) 
	{
		cout <<"Error reading Mat:  read pole "<<p<<", should be at least N_POLES = "<<maxp<<endl;
		exit(1);
	}
  
  // read mat
  if(p==maxp) fin.read( reinterpret_cast<char*>(mat), pquad*sizeof(double));
  else
	{    
		int matlen = p*p;
		streamsize readsize = psq*sizeof(double);
		streamsize skipsize = ( matlen - psq )*sizeof(double);
		
		for(int j=0; j<psq; j++)
		{
			fin.read( reinterpret_cast<char*>(mat + j*psq), readsize);
			fin.seekg( skipsize, ios::cur);
		}
  }
  
  fin.close();
  return;
	
}

//!  readMat_Formatted function 
/*! A function to read in a formatted matrix
	\param mat a matrix to read into
	\param maxp number of poles of the matrix
	\param fname a file to read the matrix from */
void readMat_Formatted(double *mat, int maxp, char* fname)
{
	
  cout << "Reading Formatted imat from "<<fname<<endl;
  const int psq = maxp * maxp;
  const int pquad = psq * psq;
	
  ifstream fin(fname);
  if(!fin)
	{
		cout <<"Error reading Mat:  cannot open "<<fname<<endl;
		exit(1);
	}
	
  // pole order
  int p;
  fin >> p ; 
  if( p < maxp) 
	{
		cout <<"Error reading Mat:  read pole "<<p<<", should be at least N_POLES = "<<maxp<<endl;
		exit(1);
	}
  
  if(p==maxp)
	{
		for(int k=0; k<pquad; k++)
			fin >> mat[k];
	}
  else
	{
		int matlen = p*p;
		double dummy;
		
		for(int j=0; j<psq; j++)
		{
			int k = j * psq;
			int i;
			for(i=0; i<psq; i++, k++)
				fin >> mat[k];
			for(i=psq; i<matlen; i++)
				fin >> dummy;
		}
	}
	
  fin.close();
  return;	
}

//!  readMats function 
/*! A function to read in an unformatted matrix
	\param iMats an array of matrices to read into
	\param p number of poles of the matrix
	\param imatpatj a path to the files to read imats from
	\param ncen the number of CG centers. Will have 1 file for each */
void readMats(vector<REAL*> &iMats, int p, char* imatpath, int ncen)
{
  iMats.resize(ncen);
  int pquad = p*p*p*p;
	
  for(int i=0; i<ncen; i++)
	{
		iMats[i] = new REAL[pquad];
		char fname[MAX_CHARNUM]; 
		int n = sprintf(fname, "%s/imat.sp%d.out.bin", imatpath ,i); 
		assert(n <= MAX_CHARNUM);
		readMat_Unformatted(iMats[i],p, fname);
	}	
}

//!  writeMat_Formatted function 
/*! A function to write out a formatted matrix
	\param mat a matrix to read into
	\param maxp number of poles of the matrix
	\param fname a file to write the matrix to */
void writeMat_Formatted(double *mat, int p, char* fname)
{
  ofstream fout(fname);
  int len = p*p;
	
  fout <<p<<endl; // pole order
  
  for(int j=0; j<len; j++)
	{
		int k = j*len;  
		for(int i=0; i<len; i++, k++)
			fout <<mat[k]<<" ";
		
		fout<<endl;
	}
	
  fout.close();
	return;
}

//!  writeExpansion function 
/*! A function to write out a multipole expansion
	\param E an expansion object to write out to file
	\param rcut a floating point of the cutoff of the expansion
	\param fname a file to write the expansion to */
void writeExpansion(const CExpan &E, REAL rcut, char* fname)
{
  ofstream fout(fname);
  fout.precision(8);
	
  fout <<E.getRange().p2()<<endl; // pole order
  fout <<E.getScale()<<endl;
  fout <<rcut<<endl;
  fout<<E<<endl;
  
  fout.close();
  return;
}

//!  readExpansion function 
/*! A subfunction to read in a multipole expansion
	\param E an expansion object to read in from file
	\param maxp number of poles of the expansion
	\param fname a file to read the expansion from 
	\param scale a floating point of the scaling factor of the expansion
	\param rcut a floating point of the cutoff of the expansion */
void readExpansion(vector<double> &V, int maxp, char* fname, REAL & scale, REAL &rcut)
{
  const int psq = maxp * maxp;
  
  V.clear();
  V.resize(psq);
  
  ifstream fin(fname);
  fin.precision(8);
	
  if(!fin)
	{
		cout <<"Error reading expansion:  cannot open "<<fname<<endl;
		exit(1);
	}
	
  // pole order
  int pRead;
  fin >> pRead ; 
  fin >> scale;
  fin >> rcut;
	
  if( pRead < maxp) 
	{
		cout <<"Warning : reading expansion:  read pole "<<pRead<<", padding zeros to N_POLES = "<<maxp<<endl;
		
		int psqRead = pRead*pRead;
		for(int k=0; k<psqRead; k++) fin >> V[k];
		for(int k=psqRead; k<psq; k++) V[k] = 0.0;
	}
  
  else
	{
		for(int k=0; k<psq; k++)  fin >> V[k];		
	}
	
  fin.close();
  return;	
}

//!  readExpansionsP function 
/*! A subfunction to read in a multipole expansion
	\param iF a multipole expansion object to read in
	\param iH a multipole expansion object to read in
	\param p number of poles of the expansion
	\param exppath a path of where the expansions are located 
	\param runname a name of the run of the expansion files 
	\param ncen an integer of the number of centers
	\param rcut a floating point of the cutoff between expansions */
void readExpansionsP(vector<CMulExpan*> &iF, vector<CMulExpan*>  &iH, 
										 int p, char* exppath, char* runname, int ncen, REAL & intraRcutoff)
{  
  iF.resize(ncen);
  iH.resize(ncen);
	
  vector<double> V;
  REAL scale;
	
  for(int ki=0; ki<ncen; ki++)
	{
		
		char fname1[MAX_CHARNUM], fname2[MAX_CHARNUM]; 
		int n;
		REAL rcut;
		
		n = sprintf(fname1, "%s/%s.%d.F.exp", exppath, runname, ki); assert(n <= MAX_CHARNUM);
		readExpansion(V, p, fname1, scale, rcut);
		if(ki==0) intraRcutoff = rcut;
		else assert( fabs(intraRcutoff-rcut) < 1e-5); 
		iF[ki] = new CMulExpan(V, CRange(p), scale);
		
		n = sprintf(fname2, "%s/%s.%d.H.exp", exppath, runname, ki); assert(n <= MAX_CHARNUM);
		readExpansion(V, p, fname2, scale, rcut);
		assert( fabs(intraRcutoff-rcut) < 1e-5); 
		iH[ki] = new CMulExpan(V, CRange(p), scale);
	}	
}

//!  generateConfigPQR function 
/*! A function to generate a PQR file from a config file
	\param configfname a character string of the config file
	\param cenfname a character string of the center file name 
	\param outfname a character string of the output file name */
void generateConfigPQR(const char* configfname, const char* cenfname, const char* outfname)
{
  vector<CPnt> rcens; 
  vector<CQuat> rots;
  readConfig(configfname, rcens, rots);
  cout <<"Read config : "<<rcens.size()<<" molecules"<<endl;
	
  // read in sphere centers and recenter them wrt molecule
  vector<CPnt> scen, posdum;
  vector<double> srad, chgdum;
  readpqr(cenfname, posdum, chgdum ,srad, scen); 
  int ncen = scen.size();
  cout <<"Read centers: "<<ncen<<" centers each"<<endl;
  CPnt initialcen;
  for(int ki=0; ki<ncen; ki++) initialcen += scen[ki];
  initialcen /= ncen;
  for(int ki=0; ki<ncen; ki++) scen[ki] -= initialcen;
	
  ofstream fout(outfname);
  for(int i=0; i<rcens.size(); i++)
	{
		CPnt rcen = rcens[i];
		CQuat rot = rots[i];
		for(int ki=0; ki < ncen; ki++)
		{
			char buf[99];
			CPnt p = (rot * scen[ki]) + rcen;
			
			
			sprintf(buf, "%6s%5d %-4s %3s %4d %8.3f%8.3f%8.3f %6.2f%6.2f\n",
							"ATOM  ",ki, "O", "CEN", ki, p.x(), p.y(),
							p.z(), 0.00, srad[ki]);
			
			fout <<buf;
			
		}//end-ki
	}//end-i
}

//!  readConfig function 
/*! A function to read in information from a config file
	\param fname a character string of the config file
	\param rcens a vector of coordinates of particle centers 
	\param rots a vector of quaternions of particle rotations */
void readConfig(const char* fname, vector<CPnt> &rcens, vector<CQuat> &rots)
{
  ifstream fin(fname);
  if(!fin)
	{
		cout <<"Error reading config:  cannot open "<<fname<<endl;
		exit(1);
	}
	
  rcens.clear();
  rots.clear();
	
  while (true)
	{
		double x, y, z, rR, rIx, rIy, rIz;
		fin >> x >> y >> z >> rR >> rIx >> rIy >> rIz;
		if(fin.eof()) break;
		rcens.push_back( CPnt(x,y,z) );
		rots.push_back( CQuat(rR, CPnt(rIx, rIy, rIz)) ); 
	}
  fin.close();	
}

////////// Unused

/*



// Read in centers and radii from centers.out file
void readCenters(const char* fname, vector<CPnt> &cens, vector<double> &radii)
{
  cens.clear();
  radii.clear();
	
  ifstream fin(fname);
  if (!fin.is_open())
	{
		cout << "Could not open file " << fname << endl;
		exit(0);
	}
	
  int j =-1;
  cout<<"Reading in centers and radius from "<<fname<<endl;
  while (!fin.eof())
	{
		REAL x,y,z,rad;
		CPnt cen;
		
		j++;
		fin>>x>>y>>z>>rad;
		cen = CPnt(x,y,z);
		cens.push_back(cen);
		radii.push_back(rad);
	}
  
  return;
}

// reads pot from APBS in openDX format
void readOpenDX(const char * fname, vector<double> &V, 
								int &nx,int &ny,int &nz, 
								double &hx,double &hy, double &hz, CPnt &pGridRef)

{
  ifstream fin(fname);
  if (!fin.is_open())
	{
		cout << "Could not open input file " << fname << endl;
		exit(0);
	}
	
  bool bParamsRead = false;
  char buf[199];
  double xmin, ymin, zmin, tempx, tempy, tempz;
  int deltaCount=0;
	
  // first get all the parameters
  while (!bParamsRead && !fin.eof())
	{
		fin.getline(buf, 199);
		
		if (strncmp(buf, "#", 1) == 0) continue;
    
		if (strncmp(buf, "object 1", 8) == 0)
		{
			sscanf(&(buf[35]), "%d %d %d", &nx, &ny, &nz);
			//	  cout <<"grid size "<<nx<<" "<<ny<<" "<<nz<<endl;
		}
		
		if (strncmp(buf, "origin", 6) == 0)
		{
			sscanf(&(buf[6]), "%le %le %le", &xmin, &ymin, &zmin);
			//  cout <<"origin " <<xmin<<" "<< ymin<<" "<< zmin<<endl;
		}
		
		if (strncmp(buf, "delta", 5) == 0)
		{
			int i;
			
			sscanf(&(buf[5]), "%le %le %le", &tempx, &tempy, &tempz);
			
			
			switch(deltaCount)
	    {
				case 0:
					hx = tempx;
					break;
				case 1:
					hy = tempy;
					break;
				case 2:
					hz = tempz;
					break;
	    }
			
			deltaCount++; 
		}
		
		if (strncmp(buf, "object 3", 8) == 0)
			bParamsRead = true;
	}
	
  //  cout <<"scale " <<hx<<" "<< hy<<" "<< hz<<endl;
	
  // then read in the potential at each point
  int n = nx*ny*nz; 
  V.resize(n);
  
  vector<double>::iterator it = V.begin(); 
  while (!fin.eof() && it != V.end() )
    fin >> *(it++);
	
  pGridRef = CPnt(xmin,ymin,zmin);
	
	
  return;
}

void 
readCompFile(char* fname, vector<CPnt> &pts, vector<REAL> &pot, int pot_column)
{
	
  pts.clear();
  pot.clear();
  
  ifstream fin(fname); 
  if (!fin.is_open())
	{
		cout << "Could not open comp file " << fname << endl;
		exit(0);
	}
	
  int id;
  REAL px, py, pz, pot1, pot2, diff;
  string s; 
  getline(fin,s);
  int ct =0;
  while (!fin.eof())
	{
		istringstream ss(s);
		ss >> id >> px >> py >> pz >> pot1 >> pot2 >> diff;
		pts.push_back( CPnt(px, py, pz));
		if(pot_column == 1) pot.push_back( pot1 );
		else if (pot_column == 2) pot.push_back( pot2 );
		else {cout <<"Error in reading compfile pot columns"<<endl; exit(1);}
		cout <<ct<<" "<<px <<" " <<py <<" "<< pz<<" "<<pot.back()<<endl;
		getline(fin, s);
		ct++;
	}
	
  fin.close();
}



void readConfigWithMolType(const char* fname, vector<CPnt> &rcens, vector<CQuat> &rots, vector<int> &moltypeid,
													 vector<int> &moltypes, vector<int> &nps)
{
  ifstream fin(fname);
  if(!fin)
	{
		cout <<"Error reading config:  cannot open "<<fname<<endl;
		exit(1);
	}
	
  rcens.clear();
  rots.clear();
  moltypeid.clear();
	
  while (true)
	{
		double x, y, z, rR, rIx, rIy, rIz;
		int  moltype;
		fin >> x >> y >> z >> rR >> rIx >> rIy >> rIz >> moltype;
		if(fin.eof()) break;
		rcens.push_back( CPnt(x,y,z) );
		rots.push_back( CQuat(rR, CPnt(rIx, rIy, rIz)) );
		moltypeid.push_back( moltype );
	}
	
  // count how many particles per moltype
  moltypes.clear();
  nps.clear();
  for (int i=0; i<moltypeid.size(); i++) 
	{
		
		bool bFound = false;
		for(int k=0; k<moltypes.size(); k++)
		{
			if(moltypeid[i] == moltypes[k])
	    {nps[k]++; bFound = true; break;}
		}
		if(!bFound)
		{
			moltypes.push_back(moltypeid[i]); 
			nps.push_back(1);
		}
	}
  
  fin.close();
	
}

*/

//!  readExpansions function 
/*! A subfunction to read in a multipole expansion
	\param iF a multipole expansion object to read in
	\param iH a multipole expansion object to read in
	\param p number of poles of the expansion
	\param exppath a path of where the expansions are located 
	\param runname a name of the run of the expansion files 
	\param ncen an integer of the number of centers
	\param rcut a floating point of the cutoff between expansions */

/*
void readExpansions(vector<CMulExpan> &iF, vector<CMulExpan>  &iH, 
										int p, char* exppath, char* runname, int ncen, REAL & intraRcutoff)
{
  iF.resize(ncen);
  iH.resize(ncen);
  
  vector<double> V;
  REAL scale;
	
  for(int ki=0; ki<ncen; ki++)
	{		
		char fname1[MAX_CHARNUM], fname2[MAX_CHARNUM]; 
		int n;
		REAL rcut;
		
		n = sprintf(fname1, "%s/%s.%d.F.exp", exppath, runname, ki); assert(n <= MAX_CHARNUM);
		readExpansion(V, p, fname1, scale, rcut);
		if(ki==0) intraRcutoff = rcut;
		else assert( fabs(intraRcutoff-rcut) < 1e-5); 
		iF[ki] = CMulExpan(V, CRange(p), scale);
		
		n = sprintf(fname2, "%s/%s.%d.H.exp", exppath, runname, ki); assert(n <= MAX_CHARNUM);
		readExpansion(V, p, fname2, scale, rcut);
		assert( fabs(intraRcutoff-rcut) < 1e-5); 
		iH[ki] = CMulExpan(V, CRange(p), scale);
	}
}

*/

