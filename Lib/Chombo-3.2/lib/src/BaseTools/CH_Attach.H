#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_ATTACH_H_
#define _CH_ATTACH_H_
#include "BaseNamespaceHeader.H"

/// call this function to register emacs-gdb to be invoked on abort.
/**
   After watching most members of the ANAG team suffer from parallel
   debugging efforts I resolved to offer something to help.

   There are two usual problems.

   1.  unless the code chokes on process rank 0, you either resort to
   printf debugging, or begin the adventure of hunting down which process
   is causing the problem and trying to use gdb 'attach' to debug the
   offending process.  One is messy and painfully slow, the other is
   finicky and difficult and you forget how to do it between times you
   need it.

   2. If you are lucky enough to actually have your code abort on process
   rank 0, you are still stuck with regular tty gdb to decipher your problem.

   All of this also depends on running your parallel process all on a
   single machine, which defeats some of the point of parallel processing.

   To address these problems, you can insert a call to 'registerDebugger()'.
   you can call it anywhere in your program.  It registers the function
   'AttachDebugger' with the ABORT signal handler.  Now, when your a process
   does something naughty and goes into abort (assert fail, MayDay, segfault,
   etc) an emacs session is launch, gdb is invoked, your binary is found (for
   symbols) and gdb attaches to your process before it has a chance to
   completely die.  The emacs window in named for the rank of the offending
   MPI process.

   Interaction with regular debug session:

   It is still perfectly fine to debug code that has called 'registerDebugger'
   in a regular gdb session, as gdb replaces signal handlers with it's own
   when it starts up your program for you.

   X11 Forwarding:

   As stated, the offending process is going to open up an emacs terminal. In
   order to do this I read the process' environment variable DISPLAY.  MPICH
   on our systems uses "ssh" to start other processes, and no amount of
   playing with mpich configure has allowed me to insert the -X command to
   enable X11 forwarding.  In addition, ssh at ANAG defaults to NOT forward
   X11. Hence, the DISPLAY environment variable for all the MPI processes
   rank>0 don't have a valid DISPLAY.  Fortunately there is an easy answer.
   Create the file ~/.ssh/config (or ~/.ssh2/config) and place the following lines in it:
      Host *
      ForwardAgent yes
      ForwardX11 yes

   This turns out to be pretty nice.  If you log into your ANAG machine from
   home using ssh -X, and then run a parallel job with mpirun you will get
   your emacs debug session forwarded from the machine the process is actually
   running on, to the machine you logged into, to your machine at home.

   If you see some message about gdb finding fd0 closed, then the failure is
   with the DISPLAY environment.  make sure you have those three lines in
   your .ssh/config

   Cavaets:

   -naturally, this approach assumes you have emacs and gdb and X11 forwarding
   on your system. If you don't then this signal handler gracefully passes the
   abort signal handler to the next handler function and you are none the
   wiser.

   -The current approach uses POSIX standard calls for it's operations,
   except: In order to find the binary for the running file I snoop into the
   /proc filesystem on linux.  This naturally will only work on linux
   operating systems.  If we find the tool useful it shouldn't be too hard to
   make it work with other OS's on a demand basis.  For those NOT running
   Linux, you are not completely without hope.  When the debugger appears, you
   will have to tell gdb where the binary is and use the 'load' command to
   fetch it and get all your debug symbols.

   -I do not know if it will be possible to have X11 forwarding working with
   batch submission systems.  It means compute nodes would need to have
   proper X client systems.  Someboday might be able to give me some pointers
   on this one.

   Another use for this feature would be babysitting very large runs.  If you
   think your big executable *might* crash then you can put this in.  If the
   code runs to completion properly, then there is not need for a debugger
   and one never gets invoked.  If your code does die, then you will find
   debugger open on your desktop attached to your program at the point where
   it failed.  Since we never seem to use core files, this might be a
   palatable option.  In parallel runs core files are just not an option.

   Tue Jan 10 13:57:20 PST 2006

   new feature:

   previously, the debugger only can be attached when the code has been
   signaled to abort, or when the mpi error handler has been called.  You get
   a debugger, attached where you want, on the MPI proc you want, but executing
   has been ended.  You cannot
   'cont'.  Additionally, if you explicitly compile the 'AttachDebugger' command
   into your code it also meant the end of the line.  no continuing.

   now, (through the POSIX, but still evil, pipe, popen, fork, etc) you can continue
   from an explicit call to AttachDebugger.  it is a two-step process.  When the
   debugger window pops up, you need to call:

       (gdb) p DebugCont()

   this should return a (void) result and put you back at the prompt. now you can
   use regular gdb commands to set break points, cont, up, fin, etc.

   You can use this feature to set an actual parallel breakpoint, you just have
   to have the patience to type your breakpoint into each window, and the fortitude
   to accept that if you make a mistake, then you have to kill all your debuggers
   and start the mpi job again.

   until I really figure out how gdbserver works, and how p4 REALLY works, this
   will have to do for now.

currently this does not work:
   In the event that the AttachDebugger cannot find a valid DISPLAY, it will
   still gdb attach to the process and read out a stacktrace to pout() so you
   can have some idea where the code died.  Unfortunately, due to a known bug
   introduced in gdb as of version 5 or so, this fallback doesn't work. There
   is a known patch for gdb, so hopefully a version of gdb will be
   available soon that works properly.  I'll still commit it as is.

   bvs
*/
int registerDebugger();  // you use this function up in main, after MPI_Init

/**
   Yeah, registerDebugger() is great if one of the Chombo CH_asserts fails or
   we call abort or exit in the code, or trip an exception, but what about when some part
   of the code calls MPI_Abort ? MPI_Abort never routes to abort but just runs to the
   exit() function after trying to take down the rest of the parallel job.  In order to
   attach to MPI_Abort calls we need to register a function with the MPI error handler
   system.  which is what this function does.  Currently registerDebugger() calls this
   function for you.
*/
int setChomboMPIErrorHandler();  // you use this function up in main, after MPI_Init

/**  Not for general consumption. you can insert this function call into
     your code and it will fire up a debugger wherever you put it.  I don't think
     the code can be continued from this point however. */
void AttachDebugger(int a_sig = 4);
#include "BaseNamespaceFooter.H"
#endif
