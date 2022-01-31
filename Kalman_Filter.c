#include "Kalman_Filter.h"
#include "fireTools.h"
#include "math.h"
#include "LQ_MPU6050_DMP.h"

#define Q_ROLL       0.0001
#define Q_PITCH      0.0001
#define R_ROLL       2.5
#define R_PITCH      2.5

extern float Pitch,Roll;
/* 说明： */
/* 1.陀螺仪y轴指向正方向，绕y轴的角是Roll角，绕x轴的角是Pitch角 */

float Pitch_Kalman,Roll_Kalman;                     //卡尔曼滤波最终输出的角度

extern signed short  aacx,aacy,aacz;	            //加速度原始数据
extern signed short  gyrox,gyroy,gyroz;             //陀螺仪原始数据
float aacx_real,aacy_real,aacz_real;                //加速度换算值
float gyrox_real,gyroy_real,gyroz_real;             //角速度换算值
float Roll_z,Pitch_z;                               //加速度计计算的角度值


float Roll_hat_pri=0,Pitch_hat_pri=0;               // 角度先验估计
float Roll_hat_pos=0,Pitch_hat_pos=0;               // 角度后验估计
float P_Roll_hat_pri=0,P_Pitch_hat_pri=0;           // 方差先验估计
float P_Roll_hat_pos=0,P_Pitch_hat_pos=0;           // 方差后验估计
float KalmanGain_Pitch,KalmanGain_Roll;             // 卡尔曼增益


float DeltaPitch=0;

float integ_angle_x=0;integ_angle_y=0;integ_angle_z=0;  //角度积分值




// ------------------------------------------------别人的代码---------------------------------------------
float Q_angle = 0.002;		//角度数据置信度，角度噪声的协方差
float Q_gyro  = 0.002;		//角速度数据置信度，角速度噪声的协方差  
float R_angle = 0.6;		//加速度计测量噪声的协方差
float dt      = 0.02;		//滤波算法计算周期，由定时器定时20ms
char  C_0     = 1;			//H矩阵值
float Q_bias, Angle_err;	//Q_bias:陀螺仪的偏差  Angle_err:角度偏量 
float PCt_0, PCt_1, E;		//计算的过程量
float K_0, K_1, t_0, t_1;	//卡尔曼增益  K_0:用于计算最优估计值  K_1:用于计算最优估计值的偏差 t_0/1:中间变量
float P[4] ={0,0,0,0};	//过程协方差矩阵的微分矩阵，中间变量
float PP[2][2] = { { 1, 0 },{ 0, 1 } };//过程协方差矩阵P
float Angle_Y_Final; 		//Y最终倾斜角度 
float Gyro_y;        		//Y轴陀螺仪数据暂存
void Kalman_Filter_Y(float Accel,float Gyro) 		
{
	Angle_Y_Final += (Gyro - Q_bias) * dt;
	P[0]=Q_angle - PP[0][1] - PP[1][0]; 
	P[1]=-PP[1][1];
	P[2]=-PP[1][1];
	P[3]=Q_gyro;	
	PP[0][0] += P[0] * dt; 
	PP[0][1] += P[1] * dt;  
	PP[1][0] += P[2] * dt;
	PP[1][1] += P[3] * dt;	
	Angle_err = Accel - Angle_Y_Final;		
	PCt_0 = C_0 * PP[0][0];
	PCt_1 = C_0 * PP[1][0];	
	E = R_angle + C_0 * PCt_0;	
	K_0 = PCt_0 / E;
	K_1 = PCt_1 / E;	
	t_0 = PCt_0;
	t_1 = C_0 * PP[0][1];
	PP[0][0] -= K_0 * t_0;		
	PP[0][1] -= K_0 * t_1;
	PP[1][0] -= K_1 * t_0;
	PP[1][1] -= K_1 * t_1;		
	Angle_Y_Final	+= K_0 * Angle_err;
	Q_bias	+= K_1 * Angle_err;	 
	Gyro_y   = Gyro - Q_bias;	 
}
//----------------------------------------------------别人的代码结束----------------------------------------



void GetAndWashData(void)
{
    // 获取原始数据
    MPU_Get_Raw_data(&aacx,&aacy,&aacz,&gyrox,&gyroy,&gyroz);
    

    /*-----------------------------------------------数据预备----------------------------------------------*/
    // 对原始数据进行单位转换,输出真实的加速度和角速度
    aacx_real = aacx/16384.0;
    aacy_real = aacy/16384.0;
    aacz_real = aacz/16384.0;
    gyrox_real = gyrox/16.4;
    gyroy_real = gyroy/16.4;
    gyroz_real = gyroz/16.4;

    // 加速度计计算(观测值)
    Pitch_z = (atan(aacy_real/(sqrt(aacx_real*aacx_real+aacz_real*aacz_real)))) * 180 / 3.1415;   // 计算y轴与水平面夹角
    Roll_z = (atan(aacx_real/(sqrt(aacy_real*aacy_real+aacz_real*aacz_real)))) * 180 / 3.1415;    // 计算x轴与水平面夹角
}

void KalmanCalculation(void)
{        
    /*-----------------------------------------卡尔曼滤波--------------------------------------------------*/
    // 参数说明：A=1，B=dt，控制量=角速度，观测值直接测得角度H=1，

    // 预测：本次角度先验估计 = 转移矩阵 * 上次角度后验估计 + 控制量
    Pitch_hat_pri = 1 * Pitch_hat_pos + gyrox_real*0.005;
    Roll_hat_pri = 1 * Roll_hat_pos +  gyroy_real*0.005;

    // 预测：本次方差先验估计 = 转移矩阵 * 上次方差后验估计 * 转移矩阵转置 + Q
    P_Pitch_hat_pri = P_Pitch_hat_pos + Q_PITCH;
    P_Roll_hat_pri = P_Roll_hat_pos + Q_ROLL;

    // 更新：计算卡尔曼增益
    KalmanGain_Pitch = (P_Pitch_hat_pri)/(P_Pitch_hat_pri+R_PITCH);
    KalmanGain_Roll = (P_Roll_hat_pri)/(P_Roll_hat_pri+R_ROLL);

    // 更新：计算后验估计
    Pitch_hat_pos = Pitch_hat_pri + KalmanGain_Pitch*(Pitch_z-Pitch_hat_pri);
    Roll_hat_pos = Roll_hat_pri + KalmanGain_Roll*(Roll_z-Roll_hat_pri);

    // 更新：计算新的方差后验估计
    P_Pitch_hat_pos = (1-KalmanGain_Pitch)*P_Pitch_hat_pri;
    P_Roll_hat_pos = (1-KalmanGain_Roll)*P_Roll_hat_pri;


    DeltaPitch = (-Pitch_hat_pos)-Pitch_Kalman;                         // 为了方便PID而进行的微分运算
    Pitch_Kalman = -Pitch_hat_pos;
    Roll_Kalman = Roll_hat_pos;

}



void TestKalman(void)
{
    //MPU6050_Init();
    GetAndWashData();
    //Kalman_Filter_Y(Pitch_z,gyroy_real);
    KalmanCalculation();
    //LQ_DMP_Init();
    //LQ_DMP_Read();
    
    // int pitchh,rolll,pitchhh,pitchzz,Angle_Y_Finall;
    // pitchh = Pitch_Kalman*100;
    // rolll = Roll_Kalman*100;
    // pitchhh = -Pitch * 100;
    // pitchzz = Pitch_z * 100;
    // Angle_Y_Finall = Angle_Y_Final * 100;
    // set_computer_value(0x02,0x01,&pitchh,1);  
    // set_computer_value(0x02,0x02,&pitchhh,1);
    // set_computer_value(0x02,0x03,&pitchzz,1);
    // set_computer_value(0x02,0x04,&Angle_Y_Finall,1);
}