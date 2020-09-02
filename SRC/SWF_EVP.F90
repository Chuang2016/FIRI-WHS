!==================================================================================================
! 子程序包: SWF_EVP (Evaporation Package for Soil Water Flow Process)
! 主要功能: 模拟地表蒸发
! 创建日期: 2020年7月3日
! 修改日志:
!          2020年8月31日, 重构了 SWF_EVP 子程序包
!==================================================================================================


SUBROUTINE SWF_EVP_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_BAS, S_IN_BAS
  USE :: SWF, ONLY : HATM
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('HATM')) THEN
    CALL RDSREA('HATM', HATM)
  ELSE
    HATM = -10000.0
  END IF

  RETURN
END SUBROUTINE SWF_EVP_DF


SUBROUTINE SWF_EVP_AL
!**************************************************************************************************
! 为数组分配内存并初始化
!**************************************************************************************************
  USE :: BAS, ONLY : NPER
  USE :: SWF, ONLY : WEPOT, WEACT
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(WEPOT)) ALLOCATE(WEPOT(NPER), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(WEACT)) ALLOCATE(WEACT(NPER), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_EVP_AL


SUBROUTINE SWF_EVP_RP
!**************************************************************************************************
! 为数组读取数据并执行准备工作
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_PER, S_IN_PER
  USE :: BAS, ONLY : NPER
  USE :: SWF, ONLY : WEPOT
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_PER 文件读取数据
  CALL RDINIT(U_IN_PER, U_LOG, S_IN_PER)
  IF (RDINQR('WEPOT')) THEN
    CALL RDFREA('WEPOT', WEPOT, NPER, NPER)
  END IF

  RETURN
END SUBROUTINE SWF_EVP_RP


SUBROUTINE SWF_EVP_FM
!**************************************************************************************************
! 计算表示地表蒸发的方程组系数
!**************************************************************************************************
  USE :: BAS, ONLY : DELV, DELT, JPER, IMAT
  USE :: SWF, ONLY : HYCONSAT, MVGA, MVGN, MVGL
  USE :: SWF, ONLY : HYCON, HNEW, RHS
  USE :: SWF, ONLY : HPONDNEW, HATM, WEPOT, QEACT
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL    :: HYCONATM, HYCON1
  REAL    :: DELV1
  REAL    :: QEMAX                ! 时段内的地表最大蒸发水量 (cm)

  I = 1
  J = IMAT(I)

  QEACT = 0.0
  CALL U_HYCON_MVG(HYCONATM, HATM, HYCONSAT(J), MVGA(J), MVGN(J), MVGL(J))
  HYCON1 = (HYCONATM + HYCON(I)) * 0.5
  DELV1 = DELV(I) * 0.5
  QEMAX = DELT * HYCON1 * ((HNEW(I) - HATM) / DELV1 - 1.0)
  QEACT = MIN(DELT * WEPOT(JPER), QEMAX)
  IF (QEACT < 0.0) QEACT = 0.0
  IF (HPONDNEW > 0.0) THEN ! 如果地表积水, 地表蒸发为水面蒸发
    HPONDNEW = HPONDNEW - QEACT
  ELSE ! 否则, 地表蒸发为土壤蒸发
    RHS(I) = RHS(I) + QEACT
  END IF

  RETURN
END SUBROUTINE SWF_EVP_FM


SUBROUTINE SWF_EVP_BD
!**************************************************************************************************
! 计算表示地表蒸发的均衡项
!**************************************************************************************************
  USE :: BAS, ONLY : JPER
  USE :: SWF, ONLY : WEACT, QEACT
  IMPLICIT NONE

  WEACT(JPER) = WEACT(JPER) + QEACT

  RETURN
END SUBROUTINE SWF_EVP_BD


SUBROUTINE SWF_EVP_DA
!**************************************************************************************************
! 释放数组内存
!**************************************************************************************************
  USE :: SWF, ONLY : WEPOT, WEACT
  IMPLICIT NONE

  IF (ALLOCATED(WEPOT))   DEALLOCATE(WEPOT)
  IF (ALLOCATED(WEACT))   DEALLOCATE(WEACT)

  RETURN
END SUBROUTINE SWF_EVP_DA

