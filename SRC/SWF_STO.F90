!==================================================================================================
! 子程序包: SWF_STO (Storage Package for Soil Water Flow Process)
! 主要功能: 模拟土壤剖面水分运动存量的变化
! 创建日期: 2020年8月30日
! 修改日志:
!==================================================================================================


SUBROUTINE SWF_STO_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE


  RETURN
END SUBROUTINE SWF_STO_DF


SUBROUTINE SWF_STO_AL
!**************************************************************************************************
! 为数组分配内存并初始化
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL
  USE :: SWF, ONLY : HINI, HOLD, HNEW, OINI, OOLD, ONEW
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(HINI)) ALLOCATE(HINI(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HOLD)) ALLOCATE(HOLD(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(HNEW)) ALLOCATE(HNEW(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OINI)) ALLOCATE(OINI(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(OOLD)) ALLOCATE(OOLD(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(ONEW)) ALLOCATE(ONEW(NCEL), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_STO_AL


SUBROUTINE SWF_STO_RP
!**************************************************************************************************
! 为数组读取数据并执行准备工作
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_SOL, S_IN_SOL
  USE :: BAS, ONLY : NCEL
  USE :: BAS, ONLY : IMAT
  USE :: SWF, ONLY : ORES, OSAT, MVGA, MVGN
  USE :: SWF, ONLY : ISWFINI
  USE :: SWF, ONLY : HINI, HOLD, HNEW, OINI, OOLD, ONEW
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 局部变量
  INTEGER :: I, J

  ! 从 S_IN_SOL 文件读取数据
  CALL RDINIT(U_IN_SOL, U_LOG, S_IN_SOL)
  IF (RDINQR('HINI')) THEN
    CALL RDFREA('HINI', HINI, NCEL, NCEL)
  ELSE
    HINI(:) = -330.0
  END IF
  IF (RDINQR('OINI')) THEN
    CALL RDFREA('OINI', OINI, NCEL, NCEL)
  ELSE
    OINI(:) = 0.2
  END IF

  ! 根据初始条件类型换算水势和含水率
  IF (ISWFINI == 0) THEN
    DO I = 1, NCEL
      J = IMAT(I)
      CALL U_O_MVG(OINI(I), HINI(I), ORES(J), OSAT(J), MVGA(J), MVGN(J))
    END DO
  ELSE
    DO I = 1, NCEL
      J = IMAT(I)
      CALL U_H_MVG(HINI(I), OINI(I), ORES(J), OSAT(J), MVGA(J), MVGN(J))
    END DO
  END IF

  HOLD(:) = HINI(:)
  HNEW(:) = HINI(:)
  OOLD(:) = OINI(:)
  ONEW(:) = OINI(:)

  RETURN
END SUBROUTINE SWF_STO_RP


SUBROUTINE SWF_STO_AD
!**************************************************************************************************
! 保存前一时段的水势和含水率
!**************************************************************************************************
  USE :: SWF, ONLY : HOLD, HNEW, OOLD, ONEW
  IMPLICIT NONE

  HOLD(:) = HNEW(:)
  OOLD(:) = ONEW(:)

  RETURN
END SUBROUTINE SWF_STO_AD


SUBROUTINE SWF_STO_FM
!**************************************************************************************************
! 计算表示土壤水分存量的方程组系数
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL
  USE :: BAS, ONLY : DELV
  USE :: SWF, ONLY : HYCAP
  USE :: SWF, ONLY : HNEW, ONEW, OOLD
  USE :: SWF, ONLY : HCOF, RHS
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I

  DO I = 1, NCEL
    HCOF(I) = HCOF(I) - DELV(I) * HYCAP(I)
    RHS(I) = RHS(I) + DELV(I) * (ONEW(I) - OOLD(I) - HYCAP(I) * HNEW(I))
  END DO

  RETURN
END SUBROUTINE SWF_STO_FM


SUBROUTINE SWF_STO_BD
!**************************************************************************************************
! 计算表示土壤水分存量的均衡项
!**************************************************************************************************
  IMPLICIT NONE


  RETURN
END SUBROUTINE SWF_STO_BD


SUBROUTINE SWF_STO_DA
!**************************************************************************************************
! 释放数组内存
!**************************************************************************************************
  USE :: SWF, ONLY : HINI, HOLD, HNEW, OINI, OOLD, ONEW
  IMPLICIT NONE

  IF (ALLOCATED(HINI)) DEALLOCATE(HINI)
  IF (ALLOCATED(HOLD)) DEALLOCATE(HOLD)
  IF (ALLOCATED(HNEW)) DEALLOCATE(HNEW)
  IF (ALLOCATED(OINI)) DEALLOCATE(OINI)
  IF (ALLOCATED(OOLD)) DEALLOCATE(OOLD)
  IF (ALLOCATED(ONEW)) DEALLOCATE(ONEW)

  RETURN
END SUBROUTINE SWF_STO_DA

