/*
       ******************************************************************
       *                                                                *
       * File:          connect_signals.c                               *
       * Author:        Edwin van der Weide                             *
       * Starting date: 03-12-2003                                      *
       * Last modified: 08-14-2005                                      *
       *                                                                *
       ******************************************************************
*/

/*
       ******************************************************************
       *                                                                *
       * This file is only compiled if signalling is supported.         *
       *                                                                *
       ******************************************************************
*/

#ifndef USE_NO_SIGNALS

#include <stdio.h>
#include <signal.h>

/*
       ******************************************************************
       *                                                                *
       * Definition of the FORTRAN functions which are called when the  *
       * user has given one of his kill signals.                        *
       *                                                                *
       ******************************************************************
*/

#ifdef FORTRAN_CAPITALS

  extern void SET_SIGNAL_WRITE(void);
  extern void SET_SIGNAL_WRITE_QUIT(void);

#elif FORTRAN_DOUBLE_UNDERSCORE

  extern void set_signal_write__(void);
  extern void set_signal_write_quit__(void);

#elif FORTRAN_NO_UNDERSCORE

  extern void set_signal_write(void);
  extern void set_signal_write_quit(void);

#else

  extern void set_signal_write_(void);
  extern void set_signal_write_quit_(void);

#endif

/*
       ******************************************************************
       *                                                                *
       * The function signal_handler is called first when the user has  *
       * given a kill signal, as this function is connected to signal,  *
       * see the function connect_signals. The function signal_handler  *
       * itself is just an interface to the FORTRAN routines.           *
       *                                                                *
       ******************************************************************
*/

void signal_handler(int sig)
{
  switch(sig)
  {
    case SIGUSR1:
#ifdef FORTRAN_CAPITALS
      SET_SIGNAL_WRITE();
#elif FORTRAN_DOUBLE_UNDERSCORE
      set_signal_write__();
#elif FORTRAN_NO_UNDERSCORE
      set_signal_write();
#else
      set_signal_write_();
#endif
      break;

    case SIGUSR2:
#ifdef FORTRAN_CAPITALS
      SET_SIGNAL_WRITE_QUIT();
#elif FORTRAN_DOUBLE_UNDERSCORE
      set_signal_write_quit__();
#elif FORTRAN_NO_UNDERSCORE
      set_signal_write_quit();
#else
      set_signal_write_quit_();
#endif
      break;
   }
}

/*
       ******************************************************************
       *                                                                *
       * The function connect_signals connects the kill commands        *
       * kill -USR1 and kill -USR2 to function signal_handler.          *
       *                                                                *
       ******************************************************************
*/

void connect_signals(void)
{
  signal(SIGUSR1, signal_handler);
  signal(SIGUSR2, signal_handler);
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void CONNECT_SIGNALS(){connect_signals();}
void connect_signals_(){connect_signals();}
void connect_signals__(){connect_signals();}

#endif  /* USE_NO_SIGNALS */
