!==================================================================================================
! 子程序包: SWF_CCF (Cell-Centered Flow Package for Soil Water Flow Process)
! 主要功能: 模拟土壤单元体之间的水分运动
! 创建日期: 2020年7月8日
! 修改日志:
!          2020年8月31日, 重构了 SWF_CCF 子程序包
!==================================================================================================


SUBROUTINE SWF_CCF_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  IMPLICIT NONE


  RETURN
END SUBROUTINE SWF_CCF_DF


SUBROUTINE SWF_CCF_AL
!**************************************************************************************************
! 为数组分配内存并初始化
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL
  USE :: SWF, ONLY : CV, QV
  IMPLICIT NONE

  IF (.NOT. ALLOCATED(CV)) ALLOCATE(CV(NCEL, 2), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(QV)) ALLOCATE(QV(NCEL, 2), SOURCE = 0.0)

  RETURN
END SUBROUTINE SWF_CCF_AL


SUBROUTINE SWF_CCF_RP
!**************************************************************************************************
! 为数组读取数据并执行准备工作
!**************************************************************************************************
  IMPLICIT NONE


  RETURN
END SUBROUTINE SWF_CCF_RP


SUBROUTINE SWF_CCF_FM
!**************************************************************************************************
! 计算表示土壤单元体之间水分运动的方程组系数
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL
  USE :: BAS, ONLY : DELV, DELT
  USE :: SWF, ONLY : HYCON, CV, RHS
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV1, DELV2
  REAL    :: HYCON1, HYCON2

  ! 第 2 至 NCEL - 1 土壤单元体
  DO I = 2, NCEL - 1
    DELV1 = (DELV(I - 1) + DELV(I)) * 0.5
    DELV2 = (DELV(I + 1) + DELV(I)) * 0.5
    HYCON1 = (DELV(I - 1) * HYCON(I - 1) + DELV(I) * HYCON(I)) / (DELV(I - 1) + DELV(I))
    HYCON2 = (DELV(I + 1) * HYCON(I + 1) + DELV(I) * HYCON(I)) / (DELV(I + 1) + DELV(I))
    CV(I, 1) = DELT * HYCON1 / DELV1
    CV(I, 2) = DELT * HYCON2 / DELV2
    RHS(I) = RHS(I) - DELT * HYCON1 + DELT * HYCON2
  END DO

  ! 第 1 土壤单元体
  I = 1
  DELV2 = (DELV(I + 1) + DELV(I)) * 0.5
  HYCON2 = (DELV(I + 1) * HYCON(I + 1) + DELV(I) * HYCON(I)) / (DELV(I + 1) + DELV(I))
  CV(I, 2) = DELT * HYCON2 / DELV2
  RHS(I) = RHS(I) + DELT * HYCON2

  ! 第 NCEL 土壤单元体
  I = NCEL
  DELV1 = (DELV(I - 1) + DELV(I)) * 0.5
  HYCON1 = (DELV(I - 1) * HYCON(I - 1) + DELV(I) * HYCON(I)) / (DELV(I - 1) + DELV(I))
  CV(I, 1) = DELT * HYCON1 / DELV1
  RHS(I) = RHS(I) - DELT * HYCON1

  RETURN
END SUBROUTINE SWF_CCF_FM


SUBROUTINE SWF_CCF_BD
!**************************************************************************************************
! 计算表示土壤单元体之间水分运动的均衡项
!**************************************************************************************************
  USE :: BAS, ONLY : NCEL
  USE :: BAS, ONLY : DELV, DELT
  USE :: SWF, ONLY : HYCON, HNEW, QV
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I
  REAL    :: DELV1, DELV2
  REAL    :: HYCON1, HYCON2

  ! 第 2 至 NCEL - 1 土壤单元体
  DO I = 2, NCEL - 1
    DELV1 = (DELV(I - 1) + DELV(I)) * 0.5
    DELV2 = (DELV(I + 1) + DELV(I)) * 0.5
    HYCON1 = (DELV(I - 1) * HYCON(I - 1) + DELV(I) * HYCON(I)) / (DELV(I - 1) + DELV(I))
    HYCON2 = (DELV(I + 1) * HYCON(I + 1) + DELV(I) * HYCON(I)) / (DELV(I + 1) + DELV(I))
    QV(I, 1) = DELT * HYCON1 * ((HNEW(I - 1) - HNEW(I)) / DELV1 + 1.0)
    QV(I, 2) = DELT * HYCON2 * ((HNEW(I + 1) - HNEW(I)) / DELV2 - 1.0)
  END DO

  ! 第 1 土壤单元体
  I = 1
  DELV2 = (DELV(I + 1) + DELV(I)) * 0.5
  HYCON2 = (DELV(I + 1) * HYCON(I + 1) + DELV(I) * HYCON(I)) / (DELV(I + 1) + DELV(I))
  QV(I, 2) = DELT * HYCON2 * ((HNEW(I + 1) - HNEW(I)) / DELV2 - 1.0)

  ! 第 NCEL 土壤单元体
  I = NCEL
  DELV1 = (DELV(I - 1) + DELV(I)) * 0.5
  HYCON1 = (DELV(I - 1) * HYCON(I - 1) + DELV(I) * HYCON(I)) / (DELV(I - 1) + DELV(I))
  QV(I, 1) = DELT * HYCON1 * ((HNEW(I - 1) - HNEW(I)) / DELV1 + 1.0)

  RETURN
END SUBROUTINE SWF_CCF_BD


SUBROUTINE SWF_CCF_DA
!**************************************************************************************************
! 释放数组内存
!**************************************************************************************************
  USE :: SWF, ONLY : CV, QV
  IMPLICIT NONE

  IF (ALLOCATED(CV)) DEALLOCATE(CV)
  IF (ALLOCATED(QV)) DEALLOCATE(QV)

  RETURN
END SUBROUTINE SWF_CCF_DA

