      INTEGER dbad_mode, dbad_phase, dbad_file, dbad_nberrors
      REAL*8 dbad_epszero, dbad_errormax, dbad_ddeps
      REAL*8 dbad_coeff, dbad_incr, dbad_sum
      REAL*8 dbad_refsum, dbad_nextrefsum
      INTEGER dbad_callsz, dbad_callindex
      PARAMETER (dbad_callsz=99)
      CHARACTER*40 dbad_callnames(dbad_callsz)
      INTEGER dbad_callcodes(dbad_callsz)
      LOGICAL dbad_calltraced(dbad_callsz)
      COMMON /DEBUG_AD_CONTROL/ dbad_mode, dbad_phase,dbad_file,
     + dbad_nberrors, dbad_epszero, dbad_errormax, dbad_ddeps,
     + dbad_coeff, dbad_incr, dbad_sum, dbad_refsum, dbad_nextrefsum,
     + dbad_callindex, dbad_callcodes, dbad_calltraced, dbad_callnames
