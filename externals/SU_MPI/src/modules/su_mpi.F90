!
!      ******************************************************************
!      *                                                                *
!      * File:          su_mpi.F90                                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-11-2003                                      *
!      * Last modified: 10-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module su_mpi
!
!      ******************************************************************
!      *                                                                *
!      * Module that contains the definition of the mpi parameters.     *
!      * Depending on the compiler flags either the module mpi or the   *
!      * include file mpif.h is loaded. In sequential mode the mpi      *
!      * functionality is faked and the mpi-parameters are simply       *
!      * defined in this module.                                        *
!      *                                                                *
!      ******************************************************************
!
#ifdef SEQUENTIAL_MODE
!
!      ******************************************************************
!      *                                                                *
!      * Part that is used when in sequential mode.                     *
!      *                                                                *
!      ******************************************************************
!
       implicit none
       save

       ! Parameter to indicate that a sequential code is built.

       logical, parameter :: SU_MPI_isSequential = .true.

       ! Parameter to indicate that parallel IO is faked. This parameter
       ! needs to be defined to avoid compiler errors. Its value does
       ! not matter.

       logical, parameter :: SU_MPI_noMPIO = .false.

       ! mpi_status

       integer, parameter :: mpi_source      = 1
       integer, parameter :: mpi_tag         = 2
       integer, parameter :: mpi_error       = 3
       integer, parameter :: mpi_status_size = 6

       ! mpi_comm

       integer, parameter :: mpi_comm_null  = 0
       integer, parameter :: mpi_comm_world = 1
       integer, parameter :: mpi_comm_self  = 2

       ! mpi_errhandler

       integer, parameter :: mpi_errhandler_null  = 0
       integer, parameter :: mpi_errors_are_fatal = 1
       integer, parameter :: mpi_errors_return    = 2

       ! mpi_group

       integer, parameter :: mpi_group_null  = 0
       integer, parameter :: mpi_group_empty = 1

       ! mpi_request

       integer, parameter :: mpi_request_null = 0

       ! mpi_op

       integer, parameter :: mpi_op_null =  0
       integer, parameter :: mpi_max     =  1
       integer, parameter :: mpi_min     =  2
       integer, parameter :: mpi_sum     =  3
       integer, parameter :: mpi_prod    =  4
       integer, parameter :: mpi_land    =  5
       integer, parameter :: mpi_band    =  6
       integer, parameter :: mpi_lor     =  7
       integer, parameter :: mpi_bor     =  8
       integer, parameter :: mpi_lxor    =  9
       integer, parameter :: mpi_bxor    = 10
       integer, parameter :: mpi_maxloc  = 11
       integer, parameter :: mpi_minloc  = 12

       ! mpi_datatype

       integer mpi_datatype_null

       integer, parameter :: mpi_char              =  1
       integer, parameter :: mpi_short             =  2
       integer, parameter :: mpi_int               =  3
       integer, parameter :: mpi_long              =  4
       integer, parameter :: mpi_unsigned_char     =  5
       integer, parameter :: mpi_unsigned_short    =  6
       integer, parameter :: mpi_unsigned          =  7
       integer, parameter :: mpi_unsigned_long     =  8
       integer, parameter :: mpi_float             =  9
       integer, parameter :: mpi_double            = 10
       integer, parameter :: mpi_long_double       = 11
       integer, parameter :: mpi_long_long         = 12
       integer, parameter :: mpi_long_long_int     = 13

       integer, parameter :: mpi_integer           = 14
       integer, parameter :: mpi_real              = 15
       integer, parameter :: mpi_double_precision  = 16
       integer, parameter :: mpi_complex           = 17
       integer, parameter :: mpi_double_complex    = 18
       integer, parameter :: mpi_logical           = 19
       integer, parameter :: mpi_character         = 20
       integer, parameter :: mpi_integer1          = 21
       integer, parameter :: mpi_integer2          = 22
       integer, parameter :: mpi_integer4          = 23
       integer, parameter :: mpi_integer8          = 24
       integer, parameter :: mpi_real4             = 25
       integer, parameter :: mpi_real8             = 26
       integer, parameter :: mpi_real16            = 27

       integer, parameter :: mpi_byte              = 28
       integer, parameter :: mpi_packed            = 29
       integer, parameter :: mpi_ub                = 30
       integer, parameter :: mpi_lb                = 31

       integer, parameter :: mpi_float_int         = 32
       integer, parameter :: mpi_double_int        = 33
       integer, parameter :: mpi_long_int          = 34
       integer, parameter :: mpi_2int              = 35
       integer, parameter :: mpi_short_int         = 36
       integer, parameter :: mpi_long_double_int   = 37

       integer, parameter :: mpi_2real             = 38
       integer, parameter :: mpi_2double_precision = 39
       integer, parameter :: mpi_2integer          = 40

       ! mpi-1 error codes and classes

       integer, parameter :: mpi_success       =  0
       integer, parameter :: mpi_err_buffer    =  1
       integer, parameter :: mpi_err_count     =  2
       integer, parameter :: mpi_err_type      =  3
       integer, parameter :: mpi_err_tag       =  4
       integer, parameter :: mpi_err_comm      =  5
       integer, parameter :: mpi_err_rank      =  6
       integer, parameter :: mpi_err_request   =  7
       integer, parameter :: mpi_err_root      =  8
       integer, parameter :: mpi_err_group     =  9
       integer, parameter :: mpi_err_op        = 10
       integer, parameter :: mpi_err_topology  = 11
       integer, parameter :: mpi_err_dims      = 12
       integer, parameter :: mpi_err_arg       = 13
       integer, parameter :: mpi_err_unknown   = 14
       integer, parameter :: mpi_err_truncate  = 15
       integer, parameter :: mpi_err_other     = 16
       integer, parameter :: mpi_err_intern    = 17
       integer, parameter :: mpi_err_in_status = 18
       integer, parameter :: mpi_err_pending   = 19

       ! mpi-2 error codes and classes

       integer, parameter :: mpi_err_access                = 28
       integer, parameter :: mpi_err_amode                 = 29
       integer, parameter :: mpi_err_assert                = 30
       integer, parameter :: mpi_err_bad_file              = 31
       integer, parameter :: mpi_err_base                  = 32
       integer, parameter :: mpi_err_conversion            = 33
       integer, parameter :: mpi_err_disp                  = 34
       integer, parameter :: mpi_err_dup_datarep           = 35
       integer, parameter :: mpi_err_file_exists           = 36
       integer, parameter :: mpi_err_file_in_use           = 37
       integer, parameter :: mpi_err_file                  = 38
       integer, parameter :: mpi_err_info_key              = 39
       integer, parameter :: mpi_err_info_nokey            = 40
       integer, parameter :: mpi_err_info_value            = 41
       integer, parameter :: mpi_err_info                  = 42
       integer, parameter :: mpi_err_io                    = 43
       integer, parameter :: mpi_err_keyval                = 44
       integer, parameter :: mpi_err_locktype              = 45
       integer, parameter :: mpi_err_name                  = 46
       integer, parameter :: mpi_err_no_mem                = 47
       integer, parameter :: mpi_err_not_same              = 48
       integer, parameter :: mpi_err_no_space              = 49
       integer, parameter :: mpi_err_no_such_file          = 50
       integer, parameter :: mpi_err_port                  = 51
       integer, parameter :: mpi_err_quota                 = 52
       integer, parameter :: mpi_err_read_only             = 53
       integer, parameter :: mpi_err_rma_conflict          = 54
       integer, parameter :: mpi_err_rma_sync              = 55
       integer, parameter :: mpi_err_service               = 56
       integer, parameter :: mpi_err_size                  = 57
       integer, parameter :: mpi_err_spawn                 = 58
       integer, parameter :: mpi_err_unsupported_datarep   = 59
       integer, parameter :: mpi_err_unsupported_operation = 60
       integer, parameter :: mpi_err_win                   = 61

       integer, parameter :: mpi_err_lastcode              = 100

       ! permanent keyvals

       integer, parameter :: mpi_keyval_invalid  = 0
       integer, parameter :: mpi_tag_ub          = 5
       integer, parameter :: mpi_host            = 6
       integer, parameter :: mpi_io              = 7
       integer, parameter :: mpi_wtime_is_global = 8

       ! results of the compare operations

       integer, parameter :: mpi_ident     = 0
       integer, parameter :: mpi_congruent = 1
       integer, parameter :: mpi_similar   = 2
       integer, parameter :: mpi_unequal   = 3

       ! topology types

       integer, parameter :: mpi_graph = 1
       integer, parameter :: mpi_cart  = 2

       ! misc constants

       integer, parameter :: mpi_max_processor_name = 256
       integer, parameter :: mpi_max_error_string   = 256
       integer, parameter :: mpi_bsend_overhead     =  32
       integer, parameter :: mpi_undefined          =  -3
       integer, parameter :: mpi_any_source         =  -2
       integer, parameter :: mpi_proc_null          =  -1
       integer, parameter :: mpi_any_tag            =  -1

       ! mpi-2

       integer, parameter :: mpi_fundamental = -1

       ! mpi-2 io

       integer, parameter :: mpi_mode_rdonly          =   2
       integer, parameter :: mpi_mode_rdwr            =   8
       integer, parameter :: mpi_mode_wronly          =   4
       integer, parameter :: mpi_mode_create          =   1
       integer, parameter :: mpi_mode_delete_on_close =  16
       integer, parameter :: mpi_mode_unique_open     =  32
       integer, parameter :: mpi_mode_excl            =  64
       integer, parameter :: mpi_mode_append          = 128
       integer, parameter :: mpi_mode_sequential      = 256

       integer, parameter :: mpi_file_null          =   0
       integer, parameter :: mpi_max_datarep_string = 128

       integer, parameter :: mpi_seek_set = 600
       integer, parameter :: mpi_seek_cur = 602
       integer, parameter :: mpi_seek_end = 604

       integer, parameter :: mpio_request_null = 0

       integer(kind=8), parameter :: mpi_displacement_current = -54278278

       integer, parameter :: mpi_order_c       = 56
       integer, parameter :: mpi_order_fortran = 57

       integer, parameter :: mpi_distribute_block     =    121
       integer, parameter :: mpi_distribute_cyclic    =    122
       integer, parameter :: mpi_distribute_none      =    123
       integer, parameter :: mpi_distribute_dflt_darg = -49767

       ! mpi-2 section 4.10

       integer, parameter :: mpi_info_null    =    0
       integer, parameter :: mpi_max_info_key =  255
       integer, parameter :: mpi_max_info_val = 1024

       ! kind values for mpi-2

       integer, parameter :: mpi_offset_kind = 8
       integer, parameter :: mpi_address_kind = 8

       ! thread-safety support levels

       integer, parameter :: mpi_thread_single     = 0
       integer, parameter :: mpi_thread_funneled   = 1
       integer, parameter :: mpi_thread_serialized = 2
       integer, parameter :: mpi_thread_multiple   = 3

       !=================================================================

       contains

         !===============================================================

         double precision function MPI_Wtime()
!
!        ****************************************************************
!        *                                                              *
!        * The subroutine MPI_Wtime returns the elapsed wall-clock time *
!        * in seconds since some time in the past. A C-function is      *
!        * called to do the actual work.                                *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         double precision :: my_time
         MPI_Wtime = my_time()

         end function MPI_Wtime

#else
!
!      ******************************************************************
!      *                                                                *
!      * Part that is used in parallel mode.                            *
!      *                                                                *
!      ******************************************************************
!
       ! Either load the include file or the module.

# ifdef USE_MPI_INCLUDE_FILE
       implicit none
       include "mpif.h"
# else
       use mpi
       implicit none
# endif

       ! Parameter to indicate that a parallel code is built.

       logical, parameter :: SU_MPI_isSequential = .false.

       ! Some definitions depending whether or not parallel IO
       ! is supported.

#  ifdef USE_NO_MPIO

       ! Parameter to indicate that parallel IO via MPIO is no
       ! supported.

       logical, parameter :: SU_MPI_noMPIO = .true.

       ! The corresponding MPI IO parameters must be defined.
       ! Otherwise SUmb will not compile.

       integer, parameter :: mpi_mode_rdonly          =   2
       integer, parameter :: mpi_mode_rdwr            =   8
       integer, parameter :: mpi_mode_wronly          =   4
       integer, parameter :: mpi_mode_create          =   1
       integer, parameter :: mpi_mode_delete_on_close =  16
       integer, parameter :: mpi_mode_unique_open     =  32
       integer, parameter :: mpi_mode_excl            =  64
       integer, parameter :: mpi_mode_append          = 128
       integer, parameter :: mpi_mode_sequential      = 256

       integer, parameter :: mpi_file_null          =   0
       integer, parameter :: mpi_max_datarep_string = 128

       integer, parameter :: mpi_seek_set = 600
       integer, parameter :: mpi_seek_cur = 602
       integer, parameter :: mpi_seek_end = 604

       integer, parameter :: mpio_request_null = 0

       integer(kind=8), parameter :: mpi_displacement_current = -54278278

       integer, parameter :: mpi_order_c       = 56
       integer, parameter :: mpi_order_fortran = 57

       integer, parameter :: mpi_distribute_block     =    121
       integer, parameter :: mpi_distribute_cyclic    =    122
       integer, parameter :: mpi_distribute_none      =    123
       integer, parameter :: mpi_distribute_dflt_darg = -49767

#  else

       ! Parameter to indicate that parallel IO via MPIO is supported.

       logical, parameter :: SU_MPI_noMPIO = .false.

#  endif

#endif

       end module su_mpi
