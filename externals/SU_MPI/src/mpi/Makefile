#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 02-21-2006                                      *
#      * Last modified: 02-24-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Makefile to create the object files of this       *
#      * directory. Note that these files should only be compiled when  *
#      * a sequential executable is created.                            *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Include the common file.                                       *
#      *                                                                *
#      ******************************************************************

SU_MPI_DIR = ../..
SU_MPI_COMMON_FILE = $(SU_MPI_DIR)/SU_MPI_Common.mk
include ${SU_MPI_COMMON_FILE}

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran and C object files to be created.         *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SU_MPI_OBJDIR)
FF90_OBJECTS  = determine_size_entity.o\
	        MPI_Abort.o\
	        MPI_Allgather.o\
	        MPI_Allgatherv.o\
	        MPI_Allreduce.o\
	        MPI_Alltoall.o\
	        MPI_Barrier.o\
	        MPI_Bcast.o\
	        MPI_Comm_rank.o\
	        MPI_Comm_size.o\
	        MPI_Comm_split.o\
	        MPI_File_close.o\
	        MPI_File_delete.o\
	        MPI_File_open.o\
	        MPI_File_read.o\
	        MPI_File_seek.o\
	        MPI_File_write.o\
	        MPI_Finalize.o\
	        MPI_Gather.o\
	        MPI_Gatherv.o\
	        MPI_Get_count.o\
	        MPI_Init.o\
	        MPI_Irecv.o\
	        MPI_Isend.o\
	        MPI_Probe.o\
	        MPI_Recv.o\
	        MPI_Reduce.o\
	        MPI_Reduce_scatter.o\
	        MPI_Send.o\
	        MPI_Waitall.o\
	        MPI_Waitany.o\
		su_mpi_terminate.o

CC_OBJECTS = MPI_File_close_c.o\
	     MPI_File_delete_c.o\
	     MPI_File_open_c.o\
	     MPI_File_read_c.o\
	     MPI_File_seek_c.o\
	     MPI_File_write_c.o\
	     MPI_Finalize_c.o\
	     MPI_Init_c.o\
	     my_time.o

OBJECTS = $(FF90_OBJECTS) $(CC_OBJECTS)

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all: 
	$(MAKE) src

clean:
	@echo "        Making clean in SU_MPI/src/mpi ... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the objects.                                      *
#      *                                                                *
#      ******************************************************************

src:	$(OBJECTS)
