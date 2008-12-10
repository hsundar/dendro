
#include "mpi.h"
#include <iostream>
#include "petsc.h"
#include "sys.h"
#include <vector>
#include "parUtils.h"
#include "testUtils.h"
#include "externVars.h"
#include "dendro.h"

class MeshRecord {
	private:
		unsigned int elementId;
		unsigned int vtxId[8];
	public:
		unsigned int getElementId() const {
			return elementId;
		}

		unsigned int getVtxId(unsigned int vtxNum) const {
			return vtxId[vtxNum];
		}

		MeshRecord() {
			this->elementId = static_cast<unsigned int>(-1);
			for(unsigned int j = 0; j < 8; j++) {
				this->vtxId[j] = static_cast<unsigned int>(-1);
			}
		}

		MeshRecord(unsigned int eId, unsigned int* vId) {
			this->elementId = eId;
			for(unsigned int j = 0; j < 8; j++) {
				this->vtxId[j] = vId[j];
			}
		}

		MeshRecord(MeshRecord const & other) {
			this->elementId = other.elementId;
			for(unsigned int j = 0; j < 8; j++) {
				((this->vtxId)[j]) = ((other.vtxId)[j]);
			}
		}

		MeshRecord & operator = (MeshRecord const  & other) {
			if( this == (& other) ) {
				return (*this);	
			}
			this->elementId = other.elementId;
			for(unsigned int j = 0; j < 8; j++) {
				((this->vtxId)[j]) = ((other.vtxId)[j]);
			}
			return (*this);
		}

		bool  operator == ( MeshRecord const  &other) const {
			bool result = (this->elementId == other.elementId);
			for(unsigned int j = 0; ( (j < 8) && result ); j++) {
				result = (result && ( ((this->vtxId)[j]) == ((other.vtxId)[j]) ) );
			}
			return result;
		}

		bool  operator < ( MeshRecord const  &other) const {
			//Sort using elementId first and then 0 thro 7 vtxId
			//Since the elements should not repeat, eId should be unique.
			if(this->elementId == other.elementId) {
				std::cout<<"WARNING: Comparing two records with the same elementId."<<std::endl;	
				unsigned int j = 0;
				while( ((this->vtxId)[j]) == ((other.vtxId)[j]) ) {
					j++;
				}
				return ( ((this->vtxId)[j]) < ((other.vtxId)[j]) );
			} else {
				return (this->elementId < other.elementId);
			}
		}

		bool  operator != (MeshRecord const  &other) const {
			return (!((*this) == other));
		}

		bool  operator > ( MeshRecord const  &other) const {
			return ( ((*this) != other) && ((*this) >= other) );
		}

		bool  operator <= ( MeshRecord const  &other) const {
			return ( ((*this) < other) || ((*this) == other) );
		}

		bool  operator >= ( MeshRecord const  &other) const {
			return (!((*this) < other));
		}

		friend std::ostream & operator << (std::ostream & os, MeshRecord const & v) ;
}; //end class MeshRecord

std::ostream & operator <<(std::ostream & os, MeshRecord const & r) {
	return (os << r.elementId <<" "<< r.vtxId[0] <<" "<<r.vtxId[1]<<" "
			<<r.vtxId[2] <<" "<< r.vtxId[3] <<" "<<r.vtxId[4]<<" "
			<<r.vtxId[5]<<" "<<r.vtxId[6]<<" "<<r.vtxId[7] );
}//end fn.

namespace par {

  template <typename T>
    class Mpi_datatype;

	template <>
		class Mpi_datatype<MeshRecord> {
			public: 
				static MPI_Datatype value() {
					static bool         first = true;
					static MPI_Datatype datatype;
					if (first) {
						first = false;
						MPI_Type_contiguous(sizeof(MeshRecord), MPI_BYTE, &datatype);
						MPI_Type_commit(&datatype);
					}
					return datatype;
				}
		};
}//end namespace

int main(int argc, char**argv) {

	PetscInitialize(&argc,&argv,0,0);
        ot::RegisterEvents();

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

	std::vector<MeshRecord> vals1;
	std::vector<MeshRecord> vals2;

	if(rank < file1Procs) {
		char fileName[100];
		sprintf(fileName,"%s_%d_%d.out",argv[1],rank,file1Procs);
		FILE * fptr = fopen(fileName,"rb");
		assert(fptr != NULL);
		unsigned int numElems = 0;
		fread(&numElems, sizeof(unsigned int), 1, fptr);
		vals1.resize(numElems);
		for(unsigned int i = 0; i < numElems; i++ ) {	
			unsigned int record[9];
			fread(record, sizeof(unsigned int), 9, fptr);
			vals1[i] = MeshRecord(record[0], (record+1));
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
			unsigned int record[9];
			fread(record, sizeof(unsigned int), 9, fptr);
			vals2[i] = MeshRecord(record[0], (record+1));
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

	bool is1uniqueAndSorted =  par::test::isUniqueAndSorted<MeshRecord>(vals1, MPI_COMM_WORLD);
	bool is2uniqueAndSorted =  par::test::isUniqueAndSorted<MeshRecord>(vals2, MPI_COMM_WORLD);

	if(!rank) {
		assert(is1uniqueAndSorted);
		assert(is2uniqueAndSorted);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"The records are sorted and unique. "<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	par::partitionW<MeshRecord>(vals1, NULL, MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Finished Partitioning List 1 "<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	par::partitionW<MeshRecord>(vals2, NULL, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Finished Partitioning List 2 "<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	{
		char outFile1[100];
		sprintf(outFile1,"MeshList1_%d_%d.out",rank,npes);
		FILE * fptr1 = fopen(outFile1,"wb");
		unsigned int numElems = (vals1.size());
		fwrite(&numElems, sizeof(unsigned int), 1, fptr1);
		for(unsigned int i = 0; i < numElems; i++ ) {	
			unsigned int record[9];
			record[0] = vals1[i].getElementId();
			for(unsigned int j = 0; j < 8; j++) {
				record[j+1] = vals1[i].getVtxId(j);
			}
			fwrite(record, sizeof(unsigned int), 9, fptr1);
		}
		fclose(fptr1);
	}

	{
		char outFile2[100];
		sprintf(outFile2,"MeshList2_%d_%d.out",rank,npes);
		FILE * fptr2 = fopen(outFile2,"wb");
		unsigned int numElems = (vals2.size());
		fwrite(&numElems, sizeof(unsigned int), 1, fptr2);
		for(unsigned int i = 0; i < numElems; i++ ) {	
			unsigned int record[9];
			record[0] = vals2[i].getElementId();
			for(unsigned int j = 0; j < 8; j++) {
				record[j+1] = vals2[i].getVtxId(j);
			}
			fwrite(record, sizeof(unsigned int), 9, fptr2);
		}
		fclose(fptr2);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank) {
		std::cout<<"Finished writing out the sorted records "<<std::endl;
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

