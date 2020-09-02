!==================================================================================================
! 农田水热盐运移数值模型 FIRI-WHS
! 机构: 中国农业科学院农田灌溉研究所
! 作者: 黄仲冬, 高阳, 张效先
! 邮箱: zdhuang@126.com
!
! 主程序 MAIN
! 创建日期: 2020年7月3日
! 修改日志:
!==================================================================================================


PROGRAM MAIN

  USE :: BAS
  USE :: SWF
  IMPLICIT NONE

  ! 局部变量
  REAL :: TSTRT, TEND

  CALL CPU_TIME(TSTRT)

  L_SWF_BAS = .FALSE.
  L_SWF_STO = .FALSE.
  L_SWF_CCF = .FALSE.
  L_SWF_GBB = .FALSE.
  L_SWF_GTB = .FALSE.
  L_SWF_INF = .FALSE.
  L_SWF_EVP = .FALSE.
  L_SWF_RWU = .FALSE.
  L_SWF_TMA = .FALSE.

  !------------------------------------------------------------------------------------------------
  ! 1. Define (DF) - 定义模型, 读取全局参数
  !------------------------------------------------------------------------------------------------
  CALL BAS_DF

  IF (L_SWF_BAS) CALL SWF_BAS_DF
  IF (L_SWF_STO) CALL SWF_STO_DF
  IF (L_SWF_CCF) CALL SWF_CCF_DF
  IF (L_SWF_GBB) CALL SWF_GBB_DF
  IF (L_SWF_GTB) CALL SWF_GTB_DF
  IF (L_SWF_INF) CALL SWF_INF_DF
  IF (L_SWF_EVP) CALL SWF_EVP_DF
  IF (L_SWF_RWU) CALL SWF_RWU_DF
  IF (L_SWF_TMA) CALL SWF_TMA_DF

  !------------------------------------------------------------------------------------------------
  ! 2. Allocate (AL) - 为数组分配内存并初始化
  !------------------------------------------------------------------------------------------------
  CALL BAS_AL

  IF (L_SWF_BAS) CALL SWF_BAS_AL
  IF (L_SWF_STO) CALL SWF_STO_AL
  IF (L_SWF_CCF) CALL SWF_CCF_AL
  IF (L_SWF_GBB) CALL SWF_GBB_AL
  IF (L_SWF_GTB) CALL SWF_GTB_AL
  IF (L_SWF_INF) CALL SWF_INF_AL
  IF (L_SWF_EVP) CALL SWF_EVP_AL
  IF (L_SWF_RWU) CALL SWF_RWU_AL
  IF (L_SWF_TMA) CALL SWF_TMA_AL

  !------------------------------------------------------------------------------------------------
  ! 3. Read & Prepare (RP) - 为数组读取数据并执行准备工作
  !------------------------------------------------------------------------------------------------
  CALL BAS_RP

  IF (L_SWF_BAS) CALL SWF_BAS_RP
  IF (L_SWF_STO) CALL SWF_STO_RP
  IF (L_SWF_CCF) CALL SWF_CCF_RP
  IF (L_SWF_GBB) CALL SWF_GBB_RP
  IF (L_SWF_GTB) CALL SWF_GTB_RP
  IF (L_SWF_INF) CALL SWF_INF_RP
  IF (L_SWF_EVP) CALL SWF_EVP_RP
  IF (L_SWF_RWU) CALL SWF_RWU_RP
  IF (L_SWF_TMA) CALL SWF_TMA_RP

  WRITE(U_LOG, "()")
  WRITE(U_LOG, "('模拟开始...')")
  WRITE(U_LOG, "('------------------------------------------------------------')")
  IF (IRUNMOD == 1) THEN
    WRITE(*, "('模拟开始...')")
    WRITE(*, "('------------------------------------------------------------')")
  END IF

  DO JPER = 1, NPER ! 应力期循环开始

    !--------------------------------------------------------------------------------------------
    ! 4. Stress (ST) - 更新随应力期变化的数据
    !--------------------------------------------------------------------------------------------
    CALL BAS_ST

    IF (L_SWF_RWU) CALL SWF_RWU_ST

    DO JSTP = 1, NSTP(JPER) ! 时段循环开始

      !--------------------------------------------------------------------------------------------
      ! 5. Advance (AD) - 更新随时段变化的数据
      !--------------------------------------------------------------------------------------------
      CALL BAS_AD

      IF (L_SWF_STO) CALL SWF_STO_AD
      IF (L_SWF_INF) CALL SWF_INF_AD

      DO JITR = 1, MXITR ! 迭代求解循环开始

        !------------------------------------------------------------------------------------------
        ! 6. Formulate (FM) - 计算土壤水分运动方程组系数
        !------------------------------------------------------------------------------------------
        IF (L_SWF_BAS) CALL SWF_BAS_FM
        IF (L_SWF_STO) CALL SWF_STO_FM
        IF (L_SWF_CCF) CALL SWF_CCF_FM
        IF (L_SWF_GBB) CALL SWF_GBB_FM
        IF (L_SWF_GTB) CALL SWF_GTB_FM
        IF (L_SWF_INF) CALL SWF_INF_FM
        IF (L_SWF_EVP) CALL SWF_EVP_FM
        IF (L_SWF_RWU) CALL SWF_RWU_FM

        !------------------------------------------------------------------------------------------
        ! 7. Approximate (AP) - 求解土壤水分运动方程组
        !------------------------------------------------------------------------------------------
        IF (L_SWF_TMA) CALL SWF_TMA_AP

        ! 如果解收敛，则退出迭代求解循环
        IF (LCONV) EXIT

      END DO ! 迭代求解循环结束

      !--------------------------------------------------------------------------------------------
      ! 8. Budget (BD) - 土壤水分运动均衡计算
      !--------------------------------------------------------------------------------------------
      IF (L_SWF_STO) CALL SWF_STO_BD
      IF (L_SWF_CCF) CALL SWF_CCF_BD
      IF (L_SWF_GBB) CALL SWF_GBB_BD
      IF (L_SWF_GTB) CALL SWF_GTB_BD
      IF (L_SWF_INF) CALL SWF_INF_BD
      IF (L_SWF_EVP) CALL SWF_EVP_BD
      IF (L_SWF_RWU) CALL SWF_RWU_BD

    END DO ! 时段循环结束

    !----------------------------------------------------------------------------------------------
    ! 9. Output (OT) - 输出模拟结果
    !----------------------------------------------------------------------------------------------
    CALL SWF_BAS_OT

    IF ((IERRCOD /= 0) .AND. (IRUNMOD /= 3)) THEN
      IF (IRUNMOD == 1) WRITE(*, "('程序停止运行!')")
      EXIT
    END IF

  END DO ! 应力期循环结束

  CALL CPU_TIME(TEND)
  WRITE(U_LOG, "('------------------------------------------------------------')")
  WRITE(U_LOG, "('模拟结束，用时 ', ES10.3, ' 秒')") TEND - TSTRT
  WRITE(U_LOG, "()")
  IF (IRUNMOD == 1) THEN
    WRITE(*, "('------------------------------------------------------------')")
    WRITE(*, "('模拟结束，用时 ', ES10.3, ' 秒')") TEND - TSTRT
    WRITE(*, *)
  END IF

  !------------------------------------------------------------------------------------------------
  ! 10. Deallocate (DA) 释放数组内存
  !------------------------------------------------------------------------------------------------
  CALL BAS_DA

  IF (L_SWF_BAS) CALL SWF_BAS_DA
  IF (L_SWF_STO) CALL SWF_STO_DA
  IF (L_SWF_CCF) CALL SWF_CCF_DA
  IF (L_SWF_GBB) CALL SWF_GBB_DA
  IF (L_SWF_GTB) CALL SWF_GTB_DA
  IF (L_SWF_INF) CALL SWF_INF_DA
  IF (L_SWF_EVP) CALL SWF_EVP_DA
  IF (L_SWF_RWU) CALL SWF_RWU_DA
  IF (L_SWF_TMA) CALL SWF_TMA_DA

  IF (IRUNMOD == 1) THEN
    WRITE(*, "('请按回车键退出...')")
    READ(*, *)
  END IF

  STOP
END PROGRAM MAIN

