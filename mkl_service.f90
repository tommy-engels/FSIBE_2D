!*******************************************************************************
!                              INTEL CONFIDENTIAL
!   Copyright(C) 1999-2008 Intel Corporation. All Rights Reserved.
!   The source code contained  or  described herein and all documents related to
!   the source code ("Material") are owned by Intel Corporation or its suppliers
!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!   suppliers and licensors. The Material contains trade secrets and proprietary
!   and  confidential  information of  Intel or its suppliers and licensors. The
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way without Intel's prior express written permission.
!   No license  under any  patent, copyright, trade secret or other intellectual
!   property right is granted to or conferred upon you by disclosure or delivery
!   of the Materials,  either expressly, by implication, inducement, estoppel or
!   otherwise.  Any  license  under  such  intellectual property  rights must be
!   express and approved by Intel in writing.
!
!*******************************************************************************
!  Content:
!      Intel(R) Math Kernel Library (MKL) FORTRAN interface for service routines
!*******************************************************************************

      INTEGER*4 MKL_ALL
      INTEGER*4 MKL_BLAS
      INTEGER*4 MKL_FFT
      INTEGER*4 MKL_VML
      INTEGER*4 MKL_DYNAMIC_TRUE
      INTEGER*4 MKL_DYNAMIC_FALSE
      PARAMETER (MKL_ALL  = 0)
      PARAMETER (MKL_BLAS = 1)
      PARAMETER (MKL_FFT  = 2)
      PARAMETER (MKL_VML  = 3)
      PARAMETER (MKL_DYNAMIC_TRUE  = 1)
      PARAMETER (MKL_DYNAMIC_FALSE = 0)

      INTERFACE
      subroutine MKLGETVERSIONSTRING(buf)
      character*(*) buf
      END
      END INTERFACE

      INTERFACE
      double precision function GETCPUFREQUENCY()
      END
      END INTERFACE

      INTERFACE
      subroutine SETCPUFREQUENCY(freq)
      double precision freq
      END
      END INTERFACE

      INTERFACE
      subroutine GETCPUCLOCKS(cpu_clocks)
      integer*8  cpu_clocks
      END
      END INTERFACE

! Threading control functions

      INTERFACE
      integer*4 function MKL_GET_MAX_THREADS()
      END
      END INTERFACE

      INTERFACE
      integer*4 function MKL_DOMAIN_GET_MAX_THREADS(domain)
      integer*4 domain
      END
      END INTERFACE

      INTERFACE
      subroutine MKL_SET_NUM_THREADS(nthrs)
      integer*4  nthrs
      END
      END INTERFACE

      INTERFACE
      integer*4 function MKL_DOMAIN_SET_NUM_THREADS(nthrs,domain)
      integer*4 nthrs
      integer*4 domain
      END
      END INTERFACE

      INTERFACE
      integer*4 function MKL_GET_DYNAMIC()
      END
      END INTERFACE

      INTERFACE
      subroutine MKL_SET_DYNAMIC(mkl_dynamic)
      integer*4 mkl_dynamic
      END
      END INTERFACE

! Memory functions

      INTERFACE
      subroutine MKL_FREE(ptr)
      pointer (ptr,mkl_dummy_var)
      integer*8 mkl_dummy_var(1)
      END
      END INTERFACE

      INTERFACE
      integer*8 function MKL_MEMSTAT(n_buff)
      integer*4 n_buff
      END
      END INTERFACE

      INTERFACE
      subroutine MKL_FREEBUFFERS()
      END
      END INTERFACE

!*******************************************************************************
