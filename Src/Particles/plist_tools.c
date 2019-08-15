/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Tools for handling the particle linked list.

  Collect various functions for adding / destroying 
  particles.         

  \authors A. Mignone (mignone@ph.unito.it)\n

  \date    March 30, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
int Particles_Insert(Particle *dummy, Data *d, char mode, Grid *grid)
/*!
 * Insert a new particle in the linked list. 
 * ...if it lies on the local processor domain.
 * Set particle id and birth rank.
 * On output return a pointer to the newly allocated particle.
 *
 *
 * \param [in]     dummy   pointer to the particle structure to be
 *                         added to the linked list.
 * \param [in,out] d       pointer to the PLUTO Data structure.
 * \param [in]     mode    specifies the action (see table).
 * \param [in]     grid    pointer to an array of Grid structures.
 *
 *  <CENTER>
 *     Action          |Ownership Check| Set Id | Insert | RestartID |
 *   ------------------|---------------|--------|------- | --------- |
 *   PARTICLES_CREATE  |  Yes          |  No    |  Yes   |   --      |
 *   PARTICLES_TRANSFER|  No           |  No    |  Yes   |   No      |
 *   PARTICLES_RESTART |  Yes          |  No    |  Yes   |   No      |
 *  </CENTER>
 *
 * Note that when mode == PARTICLES_TRANSFER we do not check the ownership
 * since a particle in a corner boundary has to be transferred twice
 * and only after the second call to Particles_Insert will actually
 * lie in the active domain.
 * 
 * \return  Return TRUE if particle was added to the list.
 *          FALSE otherwise.
 *********************************************************************** */
{
  int  nd, xBool;

/* -----------------------------------------------------------
   1. Ownership check: 
      check if the particle's coordinates fall inside current
      processor:  x(left) <= xp < x(right)  where x(left) and
      x(right) are the left and right processor interfaces 
   ----------------------------------------------------------- */
  
  xBool = TRUE;
  if (   mode == PARTICLES_CREATE
      || mode == PARTICLES_RESTART){
    DIM_LOOP(nd) {  
      if (   dummy->coord[nd] <  grid->xbeg[nd]
          || dummy->coord[nd] >= grid->xend[nd]) {
        xBool = FALSE;
      }   
    }
  }  

/* -----------------------------------------------------------
   2. If particle lies on current processor, allocate memory 
      for a new node and update linked list.
      Increment total number of particles for this task.
   ----------------------------------------------------------- */

  if (xBool){ 
    static int id = 1;
    particleNode *curr, *prev, *new_node;
    
    if (mode == PARTICLES_CREATE){
      dummy->id   = -1;     /* Set particle id     */
      dummy->tinj = g_time;     /* Set injection time to current time */
    }
    
    new_node      = malloc(sizeof(particleNode)); /* Allocate memory for a new node */
    g_usedMemory += sizeof(particleNode);
    new_node->p   = *dummy;         /* Copy input particle */
      
 /* -------------------------------------------------------
    Insert newly created node into the linked list
    (see, e.g., https://gist.github.com/mycodeschool/7429492)
    ------------------------------------------------------- */

#if 1
  /* -- Insert at head (d->PHead is always the last created node) -- */
    new_node->prev = NULL;
    if (d->PHead == NULL){   /* There're no particles. Create first one */
      d->PHead       = new_node;
      d->PHead->next = NULL;
    }else{
      new_node->next = d->PHead;
      d->PHead->prev = new_node;
      d->PHead       = new_node;  /* Change Head to new inserted particle. */
    }
#else
  /* -- Insert at tail (d->PHead is always the first created node) -- */
    if (d->PHead == NULL){   /* There're no particles. Create first one */
      d->PHead       = new_node;
      d->PHead->next = NULL;
      d->PHead->prev = NULL;
    }else{
      curr = d->PHead;
      while (curr->next != NULL) curr = curr->next;
      curr->next     = new_node;
      new_node->next = NULL;
      new_node->prev = curr;
    }
#endif
    p_nparticles++;

  }

  return xBool;
}

/* ********************************************************************* */
void Particles_Destroy(particleNode *curr, Data *d)
/*!
 *  Destroy the input node *curr and update linked list.
 *
 * \param [in]      curr  pointer to concerned particleNode that has to be
                          deleted from the list.
 * \param [in]      d     pointer to PLUTO Data structure
 
 *********************************************************************** */
{
  particleNode *prev = curr->prev;
  particleNode *next = curr->next;

  if (curr == d->PHead){  /* Remove 1st element in list */
    if (next == NULL) d->PHead = NULL;
    else{ 
      next->prev = NULL;
      d->PHead   = next;
    }  
  }else if (next == NULL){  /* Remove last element in list */
    prev->next = NULL;
  }else{
    next->prev = prev; 
    prev->next = next;
  }

  free(curr);
  g_usedMemory -= sizeof(particleNode);
  p_nparticles--;
}

/* ********************************************************************* */
void Particles_Display(Particle *p)
/*!
 * Prints the relevant quantities for particular particle.
 *
 * \param [in]      p  pointer to Particle structure whose values are to be 
                       printed.
 *********************************************************************** */
{
  print ("////////////////////////////////////////////////////////////////\n");
  print ("// id:           %d\n", p->id);
  print ("// (x1,x2,x3):   %f, %f, %f\n", p->coord[IDIR],  p->coord[JDIR], p->coord[KDIR]);
  print ("// (v1,v2,v3):   %f, %f, %f\n", p->speed[IDIR], p->speed[JDIR], p->speed[KDIR]);
  print ("// (i,j,k):      %d, %d, %d\n", p->cell[IDIR],  p->cell[JDIR], p->cell[KDIR]);
  print ("// tinj:         %f\n", p->tinj);
  print ("////////////////////////////////////////////////////////////////\n"); 
}

/* ********************************************************************* */
int Particles_Number(particleNode* PNodeRef)
/*! 
 *  Count and return the number of particles owned by the current 
 *  process rnk.
 *
 * \param [in]  PNodeRef  pointer to headnode of the list 
                          whose count has to be estimated.
 *
 *********************************************************************** */
{
  particleNode* CurNode = PNodeRef;

  int cnt = 0;
  while(CurNode != NULL){
    cnt++;
    CurNode = CurNode->next;
  }
  return cnt;
}

/* ********************************************************************* */
void Particles_ShowList(particleNode *PHead, int show_particles)
/*! 
 *  Prints the memory addresses (and relevant quantities of particles) 
 *  of the linked list.
 *
 *  \param [in]  PHead              pointer to headnode of the particle list.
 *  \param [in]  show_particles     An integer value, if set to 1 this routine
 *                                  will also display relevant quantities for
 *                                  every particle in the list.
 *********************************************************************** */
{
  int count = 0;
  particleNode *curr = PHead;
  
  while (curr != NULL){
    print ("# %02d, curr = %p  (prev = %p, next = %p)\n", 
           count, curr, curr->prev, curr->next);
    if (show_particles) Particles_Display(&(curr->p));       
    count++;   
    curr = curr->next;
  }
  print ("------- Done ----- \n");
}
