
#include "mpi.h"
#include <iostream>
#include "petsc.h"
#include "sys.h"
#include <vector>
#include "parUtils.h"
#include "externVars.h"
#include "dendro.h"

class Vertex {
	private:
		unsigned int x, y, z, id;
	public:
		unsigned int getX() const {
			return x;
		}

		unsigned int getY() const {
			return y;
		}

		unsigned int getZ() const {
			return z;
		}

		unsigned int getId() const {
			return id;
		}

		Vertex() {
			this->x = static_cast<unsigned int>(-1);
			this->y = static_cast<unsigned int>(-1);
			this->z = static_cast<unsigned int>(-1);
			this->id = static_cast<unsigned int>(-1);
		}

		Vertex(unsigned int xVal, unsigned int yVal, unsigned int zVal, unsigned int idVal) {
			this->x = xVal;
			this->y = yVal;
			this->z = zVal;
			this->id = idVal;
		}

		Vertex(Vertex const & other) {
			this->x = other.x;
			this->y = other.y;
			this->z = other.z;
			this->id = other.id;
		}

		Vertex & operator = (Vertex const  & other) {
			if( this == (&other) ) {
				return (*this);	
			}
			this->x = other.x;
			this->y = other.y;
			this->z = other.z;
			this->id = other.id;
			return (*this);
		}

		bool  operator == ( Vertex const  &other) const {
			return( (this->x == other.x) && (this->y == other.y) && 
					(this->z == other.z) && (this->id == other.id) );
		}

		bool  operator != (Vertex const  &other) const {
			return (!((*this) == other));
		}

		bool  operator < ( Vertex const  &other) const {
			//Sort using x first and then  y and then z and then id
			if(this->x == other.x) {
				if( this->y == other.y) {
					if( this->z == other.z) {
						return (this->id < other.id);
					}else {
						return (this->z < other.z);
					}
				}else {
					return(this->y < other.y);
				}
			} else {
				return (this->x < other.x);
			}
		}

		bool  operator > ( Vertex const  &other) const {
			return ( ((*this) != other) && ((*this) >= other) );
		}

		bool  operator <= ( Vertex const  &other) const {
			return ( ((*this) < other) || ((*this) == other) );
		}

		bool  operator >= ( Vertex const  &other) const {
			return (!((*this) < other));
		}

		friend std::ostream & operator << (std::ostream & os, Vertex const & v) ;
}; //end class Vertex

std::ostream & operator <<(std::ostream & os, Vertex const & v) {
	return (os << v.x <<" "<< v.y <<" "<<v.z<<" "<<v.id);
}//end fn.

namespace par {

  template <typename T>
    class Mpi_datatype;
  
	template <>
		class Mpi_datatype<Vertex> {
			public: 
				static MPI_Datatype value() {
					static bool         first = true;
					static MPI_Datatype datatype;
					if (first) {
						first = false;
						MPI_Type_contiguous(sizeof(Vertex), MPI_BYTE, &datatype);
						MPI_Type_commit(&datatype);
					}
					return datatype;
				}
		};
}//end namespace

int main(int argc, char**argv) {

	PetscInitialize(&argc,&argv,0,0);
        ot::RegisterEvents();//Register OTK Functions

	if(argc < 4) {
		std::cout<<"exe fileBase file1Procs file2Procs "<<std::endl;
		assert(false);
	}

	int rank, npes;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);

	unsigned int file1Procs = atoi(argv[2]);
	unsigned int file2Procs = atoi(argv[3]);

	assert(file1Procs <= file2Procs);
	assert(file2Procs <= npes);

	std::vector<Vertex> vals1;
	std::vector<Vertex> vals2;

	if(rank < file1Procs) {
		char fileName[100];
		sprintf(fileName,"%s_%d_%d.out",argv[1],rank,file1Procs);
		FILE * fptr = fopen(fileName,"rb");
		assert(fptr != NULL);
		unsigned int numElems = 0;
		fread(&numElems, sizeof(unsigned int), 1, fptr);
		vals1.resize(numElems);
		for(unsigned int i = 0; i < numElems; i++ ) {	
			unsigned int xyzi[4];
			fread(xyzi, sizeof(unsigned int), 4, fptr);
			vals1[i] = Vertex(xyzi[0], xyzi[1], xyzi[2], xyzi[3]);
		}
		fclose(fptr);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Read File1 Successfully."<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank < file2Procs) {
		char fileName[100];
		sprintf(fileName,"%s_%d_%d.out",argv[1],rank,file2Procs);
		FILE * fptr = fopen(fileName,"rb");
		assert(fptr != NULL);
		unsigned int numElems = 0;
		fread(&numElems, sizeof(unsigned int), 1, fptr);
		vals2.resize(numElems);
		for(unsigned int i = 0; i < numElems; i++ ) {	
			unsigned int xyzi[4];
			fread(xyzi, sizeof(unsigned int), 4, fptr);
			vals2[i] = Vertex(xyzi[0], xyzi[1], xyzi[2], xyzi[3]);
		}
		fclose(fptr);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Read File2 Successfully."<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	DendroIntL locSz1 = vals1.size();
	DendroIntL globSz1 = 0;
	par::Mpi_Reduce<DendroIntL>(&locSz1, &globSz1, 1, MPI_SUM, 0, MPI_COMM_WORLD);

	DendroIntL locSz2 = vals2.size();
	DendroIntL globSz2 = 0;
	par::Mpi_Reduce<DendroIntL>(&locSz2, &globSz2, 1, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		assert(globSz1 == globSz2);
		std::cout << "Global Size: " << globSz1 << std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<Vertex> vals1Tmp;
	std::vector<Vertex> vals2Tmp;

	par::sampleSort<Vertex>(vals1, vals1Tmp, MPI_COMM_WORLD); 
	par::sampleSort<Vertex>(vals2, vals2Tmp, MPI_COMM_WORLD); 

	vals1 = vals1Tmp;
	vals2 = vals2Tmp;
	vals1Tmp.clear();
	vals2Tmp.clear();

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Finished Sorting Vertices "<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	par::partitionW<Vertex>(vals1, NULL, MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Finished Partitioning List 1 "<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	par::partitionW<Vertex>(vals2, NULL, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Finished Partitioning List 2 "<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	{
		char outFile1[100];
		sprintf(outFile1,"VtxList1_%d_%d.out",rank,npes);
		FILE * fptr1 = fopen(outFile1,"wb");
		unsigned int numElems = (vals1.size());
		fwrite(&numElems, sizeof(unsigned int), 1, fptr1);
		for(unsigned int i = 0; i < numElems; i++ ) {	
			unsigned int record[4];
			record[0] = vals1[i].getX();
			record[1] = vals1[i].getY();
			record[2] = vals1[i].getZ();
			record[3] = vals1[i].getId();
			fwrite(record, sizeof(unsigned int), 4, fptr1);
		}
		fclose(fptr1);
	}

	{
		char outFile2[100];
		sprintf(outFile2,"VtxList2_%d_%d.out",rank,npes);
		FILE * fptr2 = fopen(outFile2,"wb");
		unsigned int numElems = (vals2.size());
		fwrite(&numElems, sizeof(unsigned int), 1, fptr2);
		for(unsigned int i = 0; i < numElems; i++ ) {	
			unsigned int record[4];
			record[0] = vals2[i].getX();
			record[1] = vals2[i].getY();
			record[2] = vals2[i].getZ();
			record[3] = vals2[i].getId();
			fwrite(record, sizeof(unsigned int), 4, fptr2);
		}
		fclose(fptr2);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Finished writing out the sorted vertices "<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	assert(vals1.size() == vals2.size());

	MPI_Barrier(MPI_COMM_WORLD);

	bool Ipassed = false;
	if(vals1 == vals2) {
		Ipassed = true;
	}

	bool *IpassedList = new bool[npes];
        par::Mpi_Allgather<bool>(&Ipassed, IpassedList, 1, MPI_COMM_WORLD);

	bool allPassed = true;
	unsigned int lastP = 0;
	for(unsigned int i = 0; ((i < npes) && allPassed); i++) {
		allPassed = (allPassed && IpassedList[i]);
		lastP = i;
	}
	delete [] IpassedList;

	if(!allPassed) {
		if(rank == lastP) {
			std::cout<<"First Failing on Proc: "<< lastP <<std::endl;
			for(unsigned int i = 0; i < vals1.size(); i++) {
				if(vals1[i] != vals2[i]) {
					std::cout<<"i = "<<i<<" val1 = "<<vals1[i]<<" val2 = "<<vals2[i]<<std::endl;
					break;	
				}
			}
		}
	}

	vals1.clear();
	vals2.clear();

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		if(allPassed) {
			std::cout<<"The files match."<<std::endl;
		}else {
			std::cout<<"The files do not match."<<std::endl;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	PetscFinalize();

}

