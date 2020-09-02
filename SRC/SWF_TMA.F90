!==================================================================================================
! 子程序包: SWF_TMA (Tridiagonal Matrix Algorithm Package for Soil Water Flow Process)
! 主要功能: 采用 Thomas 算法求解土壤水分运动有限差分方程组
! 创建日期: 2020年7月7日
! 修改日志:
!          2020年8月31日, 重构了 SWF_TMA 子程序包
!==================================================================================================


SUBROUTINE SWF_TMA_DF
!**************************************************************************************************
! 读取模型定义参数
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG, U_IN_BAS, S_IN_BAS
  USE :: BAS, ONLY : MXITR
  USE :: SWF, ONLY : HTOL, OTOL
  IMPLICIT NONE

  ! 函数声明
  LOGICAL :: RDINQR

  ! 从 S_IN_BAS 文件读取数据
  CALL RDINIT(U_IN_BAS, U_LOG, S_IN_BAS)
  IF (RDINQR('MXITR')) THEN
    CALL RDSINT('MXITR', MXITR)
  ELSE
    MXITR = 30
  END IF
  IF (RDINQR('HTOL')) THEN
    CALL RDSREA('HTOL', HTOL)
  ELSE
    HTOL = 1.0
  END IF
  IF (RDINQR('OTOL')) THEN
    CALL RDSREA('OTOL', OTOL)
  ELSE
    OTOL = 0.0001
  END IF

  RETURN
END SUBROUTINE SWF_TMA_DF


SUBROUTINE SWF_TMA_AL
!**************************************************************************************************
! 为数组分配内存并初始化
!**************************************************************************************************
  IMPLICIT NONE


  RETURN
END SUBROUTINE SWF_TMA_AL


SUBROUTINE SWF_TMA_RP
!**************************************************************************************************
! 为数组读取数据并执行准备工作
!**************************************************************************************************
  IMPLICIT NONE


  RETURN
END SUBROUTINE SWF_TMA_RP


SUBROUTINE SWF_TMA_AP
!**************************************************************************************************
! 求解土壤水分运动有限差分方程组
!**************************************************************************************************
  USE :: BAS, ONLY : U_LOG
  USE :: BAS, ONLY : IRUNMOD, IERRCOD
  USE :: BAS, ONLY : NCEL, IMAT
  USE :: BAS, ONLY : MXITR, JITR, LCONV
  USE :: SWF, ONLY : ORES, OSAT
  USE :: SWF, ONLY : MVGA, MVGN
  USE :: SWF, ONLY : HNEW, ONEW
  USE :: SWF, ONLY : HTOL, OTOL
  USE :: SWF, ONLY : CV, HCOF, RHS
  IMPLICIT NONE

  ! 局部变量
  INTEGER :: I, J
  REAL :: BET
  REAL, DIMENSION(:), ALLOCATABLE :: D
  REAL, DIMENSION(:), ALLOCATABLE :: E
  REAL, DIMENSION(:), ALLOCATABLE :: F
  REAL, DIMENSION(:), ALLOCATABLE :: R
  REAL, DIMENSION(:), ALLOCATABLE :: G
  REAL, DIMENSION(:), ALLOCATABLE :: U
  REAL, DIMENSION(:), ALLOCATABLE :: V

  ! 为局部变量分配内存
  IF (.NOT. ALLOCATED(D)) ALLOCATE(D(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(E)) ALLOCATE(E(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(F)) ALLOCATE(F(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(R)) ALLOCATE(R(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(G)) ALLOCATE(G(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(U)) ALLOCATE(U(NCEL), SOURCE = 0.0)
  IF (.NOT. ALLOCATED(V)) ALLOCATE(V(NCEL), SOURCE = 0.0)

  ! 计算三对角线性方程组的系数和右端项
  D(:) = CV(:, 1)
  E(:) = - CV(:, 1) - CV(:, 2) + HCOF(:)
  F(:) = CV(:, 2)
  R(:) = RHS(:)

  ! 求解三对角线性方程组
  BET = E(1)
  IF (ABS(BET) < TINY(0.0)) THEN
    WRITE(U_LOG, "('I = ', I12,' BETA = ', ES12.3)") 1, BET
    WRITE(U_LOG, "('错误 1 -- 土壤水分运动方程组系数错误!')")
    IF (IRUNMOD == 1) WRITE(*, "('I = ', I12,' BETA = ', ES12.3)") 1, BET
    IF (IRUNMOD == 1) WRITE(*, "('错误 1 -- 土壤水分运动方程组系数错误!')")
    IERRCOD = 1
  END IF
  U(1) = R(1) / BET
  DO I = 2, NCEL
    G(I) = F(I - 1) / BET
    BET = E(I) - D(I) * G(I)
    IF (ABS(BET) < TINY(0.0)) THEN
      WRITE(U_LOG, "('I = ', I12,' BETA = ', ES12.3)") I, BET
      WRITE(U_LOG, "('错误 2 -- 土壤水分运动方程组系数错误!')")
      IF (IRUNMOD == 1) WRITE(*, "('I = ', I12,' BETA = ', ES12.3)") I, BET
      IF (IRUNMOD == 1) WRITE(*, "('错误 2 -- 土壤水分运动方程组系数错误!')")
      IERRCOD = 1
    END IF
    U(I) = (R(I) - D(I) * U(I - 1)) / BET
  END DO
  DO I = NCEL - 1, 1, -1
    U(I) = U(I) - G(I + 1) * U(I + 1)
  END DO

  DO I = 1, NCEL
    J = IMAT(I)
    CALL U_O_MVG(V(I), U(I), ORES(J), OSAT(J), MVGA(J), MVGN(J))
  END DO

  ! 判断解是否收敛
  LCONV = .FALSE.
  IF ((MAXVAL(ABS(HNEW(:) - U(:))) < HTOL) .AND. (MAXVAL(ABS(ONEW(:) - V(:))) < OTOL)) THEN
    LCONV = .TRUE.
  END IF

  ! 如果不收敛, 输出一个警告信息
  IF ((JITR >= MXITR) .AND. (.NOT. LCONV)) THEN
    WRITE(U_LOG, "('警告 -- 土壤水分运动方程组的解不收敛!')")
    IF (IRUNMOD == 1) WRITE(*, "('警告 -- 土壤水分运动方程组的解不收敛!')")
    IERRCOD = 0
  END IF

  ! 保存方程组的解
  HNEW(:) = U(:)
  ONEW(:) = V(:)

  ! 释放局部变量的内存
  IF (ALLOCATED(D)) DEALLOCATE(D)
  IF (ALLOCATED(E)) DEALLOCATE(E)
  IF (ALLOCATED(F)) DEALLOCATE(F)
  IF (ALLOCATED(R)) DEALLOCATE(R)
  IF (ALLOCATED(G)) DEALLOCATE(G)
  IF (ALLOCATED(U)) DEALLOCATE(U)
  IF (ALLOCATED(V)) DEALLOCATE(V)

  RETURN
END SUBROUTINE SWF_TMA_AP


SUBROUTINE SWF_TMA_DA
!**************************************************************************************************
! 释放数组内存
!**************************************************************************************************
  IMPLICIT NONE

  RETURN
END SUBROUTINE SWF_TMA_DA