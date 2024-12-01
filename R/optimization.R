# \arguments{
#   \item{S0}{An array recording the initial prices of different assets(eg:S0 <- c(20804.11,5969.34,2406.67))}
#   \item{N}{Number of Brownian motion predictions}
#   \item{ST}{The end-of-period price of each asset is predicted N times by Brownian motion, with num_w rows and N columns}
#   \item{alpha}{alpha}
#   \item{num_w}{Number of assets}
#   \item{num_p}{The number of options on each asset, an array of length num_w}
#   \item{C}{The current market price of the option, with num_w rows and sum_p columns}
#   \item{K}{Strike price for each asset, with num_w rows and sum_p columns}
#   \item{rf}{Risk-free rate of return}
#   \item{T}{Duration of asset price change research, using days/365}
#   \item{Mu}{Each asset Mu in the Brownian motion, an array of length num_w}
#   \item{Sigm}{The volatility of each asset in the Brownian motion, an array of length num_w}
# }
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'



#Other:
#The data is in data_df.rda and data_k.rda. data_df corresponds to parameter C, which is the current market price of the option.
#data_K corresponds to parameter K, which is the option's expiration exercise price
#C<- as.matrix(data_df)
#K<- as.matrix(data_k)
#An Example is S0 <- c(20804.11,5969.34,2406.67), alpha <- 0.95, N <- 100, num_w <- 3, num_p <- c(22,95,51)

#C<- as.matrix(data_df)
#K<- as.matrix(data_k)
optimization<-function(S0,N,alpha,num_w,num_p,C,K,rf,T,Mu,Sigm) {
  p1 <- num_p[1]
  p2 <- num_p[2]
  p3 <- num_p[3]
  library(lpSolve)

  #运用几何布朗运动 生成ST
  set.seed(123)
  generate_gbm_st <- function(s0, mu, sigm, T, n = 1) {
    # 参数说明：
    # S0    - 初始价格
    # mu    - 漂移率
    # sigma - 波动率
    # T     - 总时间
    # n     - 模拟路径数量（默认 1）
    # 生成 n 个标准正态随机数
    Z <- rnorm(n, mean = 0, sd = 1)

    # 计算 S_T
    ST0 <- s0 * exp((mu - 0.5 * sigm^2) * T + sigm * sqrt(T) * Z)

    # 返回结果
    return(ST0)
  }
  ST <- matrix(0,nrow=num_w,ncol=N)

  for(j in (1:num_w)){
    for( i in (1:N)){
      ST[j,i] <- generate_gbm_st(s0=S0[j],mu= Mu[j] ,sigm= Sigm[j] ,T= 7/365,n=1)
    }
  }

  ###线性优化
  ##目标系数
  f.obj <- c(1,rep(1/(N*(1-alpha)),N),rep(0,num_w),rep(0,sum_p))
  total_cols <- length(f.obj)

  ###约束矩阵定义（列）
  # A[1] : q
  # A[2:N+1] : z
  # A[N+2:N+1+num_w] : w
  # A[N+2+num_w:N+1+num_w+p1] :p1l
  # A[N+2+num_w+p1:N+1+num_w+p1+p2] :p2l
  # A[N+2+num_w+p1+p2:N+1+num_w+p1+p2+p3] :p3l

  ###约束矩阵定义（行）
  # A[1:num_w+1]: 约束 1
  # A[num_w+2:num_w+1+sum_p]:约束 2
  # ...

  ## 约束1：wj ≥ 0, ∑wj = 1
  A1 <- matrix(0, nrow = num_w+1 , ncol = total_cols)

  for (i in 1:(num_w)) {
    A1[i, (i + N+2-1)] <- 1
    A1[(num_w+1),(i+N+2-1)] <-1
  }
  B1 <- c(rep(0, num_w), 1)
  direction1 <- c(rep(">=", num_w), "==")

  ## 约束2：pjl >= 0
  A2 <- matrix(0, nrow = sum_p , ncol = total_cols)
  for(i in 1:(sum_p)){
    A2[i, (N+2+num_w-1+ i)] <-1
  }
  B2 <- c(rep(0,sum_p))
  direction2 <- c(rep(">=",sum_p))

  ##() 约束 3：对于每一个j，pjl的和 <= (wj/S0j) >>>> -wj+S0*pjl的和<=0
  A3 <- matrix(0,nrow = num_w, ncol =total_cols)
  ##因为只有三个资产，直接展开写了
  for(i in (N+2+num_w):(N+1+num_w+p1)){
    A3[1,i] <- 1
  }
  for(i in (N+2+num_w+p1):(N+1+num_w+p1+p2)){
    A3[2,i] <- 1
  }
  for(i in (N+2+num_w+p1+p2):(N+1+num_w+p1+p2+p3)){
    A3[3,i] <- 1
  }
  for (i in 1:(num_w)) {
    A3[i, (i + N+2-1)] <- -1/S0[i]
  }
  B3 <- c(rep(0,num_w))
  direction3 <- c(rep("<=",num_w))

  ## 约束 4：
  ##计算中间量
  tem <- array(0,dim=c(num_w,(max(num_p)),N))
  tem_i <- matrix(0,nrow=num_w,ncol=(max(num_p)))  ##对不同 i 求和
  for(j in (1:num_w)){
    for(l in (1:num_p[j])){
      for(i in (1:N)){
        tem[j,l,i] <- C[j,l] * exp(rf * T)-max((ST[j,i]-K[j,l]),0)
        #cat("\nj,l,i,max((ST[j,i]-K[j,l]),0)",j,l,i,ST[j,i]-K[j,l])
        tem_i[j,l] <- tem_i[j,l]+tem[j,l,i]
        #print(tem_i[j,l])
      }
    }
  }
  ###填p的系数
  A4 <- matrix(0,nrow = 1, ncol =total_cols)
  count_tail <- N+2+num_w
  for(j in (1:num_w)){
    for(l in (1:num_p[j])){
      A4[1,(count_tail)] <- tem_i[j,l]
      count_tail <- count_tail+1
    }
  }
  ###对 ST在i上加和
  ST_i <- rep(0,num_w)
  for(j in (1:num_w)){
    for(i in (1:N)){
      ST_i[j] <- ST_i[j]+ST[j,i]
    }
  }
  ##填w的系数
  count_tail2 <- N+2
  for(j in(1:(num_w))){
    A4[1,count_tail2] <- ST_i[j]/S0[j]
    count_tail2 <- count_tail2 + 1
  }
  r_target <- 0.007
  B4 <- c(N*(r_target+1))
  direction4 <- c(">=")

  ## 约束 5：-ri-zi-q <= 0 >>>>
  A5 <- matrix(0,nrow = N, ncol = total_cols)
  for(i in (1:N)){
    A5[i,1] <- 1   ## q
    A5[i,i+1] <- 1  ## z
    count_tail3 <- N+2
    count_tail4 <- N+2+num_w
    for(j in (1:num_w)){
      A5[i,count_tail3] <- ST[j,i]/S0[j] ##填w系数
      count_tail3 <- count_tail3 + 1
      for(l in (1:num_p[j])){
        A5[i,count_tail4] <- tem[j,l,i] ##填p系数
        count_tail4 <- count_tail4 + 1
      }
    }
  }
  B5 <- c(rep(1,N))
  direction5 <- c(rep(">=",N))
  ## 约束 6：zi >= 0
  A6 <- matrix(0,nrow = N, ncol =total_cols)
  for(i in (1:N)){
    A6[i,(2-1+i)] <-1
  }
  B6 <-rep(0,N)
  direction6 <- c(rep(">=",N))

  ##合并
  A <- rbind(A1, A2, A3, A4, A5, A6)
  #f.con <- Ax
  Direction <- c(direction1,direction2,direction3,direction4,direction5,direction6)
  rt <- seq(0.003,0.009,by = 0.00008)
  ##输出
  answer <- c(rep(0,272))
  riskanswer <- c(0)
  returnanswer <- c(0)
  for(rtarget1 in rt){
    B4 <- c(N*(rtarget1+1))
    #print(B4)
    B <- c(B1, B2, B3, B4, B5, B6)
    result <- lp("min", f.obj, const.mat = A, const.rhs = B, const.dir = Direction)
    # print(result[["solution"]][102])
    # print(result[["solution"]][103])
    # print(result[["solution"]][104])
    # cat("##############################\n")

    ##保留每次优化的结果
    answer<- rbind(answer,result[["solution"]])
    riskanswer <- c(riskanswer,result$objval)

    ##计算 Return
    solution <- result[["solution"]]
    RI <- matrix(0,nrow = N, ncol = total_cols)
    for(i in (1:N)){
      count_tail3 <- N+2
      count_tail4 <- N+2+num_w
      for(j in (1:num_w)){
        RI[i,count_tail3] <- ST[j,i]/S0[j] ##填w系数
        count_tail3 <- count_tail3 + 1
        for(l in (1:num_p[j])){
          RI[i,count_tail4] <- tem[j,l,i] ##填p系数
          count_tail4 <- count_tail4 + 1
        }
      }
    }
    test1 <- RI * solution
    Return <- (sum(test1)-N)/N
    returnanswer <- c(returnanswer,Return)
  }
  return(result)
}
