c********1*********1*********1*********1*********1*********1********1**
c Pour ETACHA4
c.....calcula la secc. eficaz total de excitacion desde los
c.....estados ns, np1 y np-1 al estado final nlm, por S.Eikonal
c.....Autor: Cesar A. Ramirez            Fecha: 10-12-2012
c         zt=carga del nucleo blanco
c         zp=carga del nucleo proyectil 
c         ei=energia del estado ligado inicial
c         ef=energia del estado ligado final
c         n,l,m = numeros cuanticos del estado final
c         eta=componente normal del momento transferido
c**********************************************************************
c----- Program to calculte total cross sections of  
c----- monoelectronic atoms excitation, for colision with nude projectiles 
c----- in the SE (Symmetric-Eikonal) aproximation.
c--------------------------------------------------
c......changes in the present version (JPR 02/2012)
c----- zp: Projectile charge (to be excited)
c----- zt: Target charge
c......SC and ASC corrections
c-------------------------------------------------
c----- (n0,l0,m0) quantum numbers of initial state 
c----- (n,l,m) quantum numbers of final state
c----- SUBROUTINES:
c----- * rns, rnp0, rnp1, rnp1m: excitation from ns, np0, np1 np-1 to n'lm
c-----                    respectively with n'lm arbitraries (n'>n)
c----- * intdef: for definite integrals from 0 to infinity
c----- * fac: factorial function
c----- * cg: Clebsch-Gordan coeficients
c----- * armonic: spherical harmonic function
c----- * plgndr: Legendre polinomials
c----- * gamac: gama function
c----- * cf21d, chyper, hypgfx: hypergeometric functions F21(ca,cb,cc,cx)
c**********************************************************************
      subroutine SEnlm(zp8,zt8,E8)
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      real*4 SeSE,StSE
	character kk*2
      end if
      if(dabs(csum).gt.preci*dabs(cint))then
      go to 39
      else
      continue
      end if            
      return
      end
c***********************************************************************