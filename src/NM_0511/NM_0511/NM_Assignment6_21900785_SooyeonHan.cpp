/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Sooyeon Han
Created          : 18-05-2021
Modified         : 18-05-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	6		// enter your assignment number
#define eval		0		// set 0

#include "../../../include/myNM.h"

# include "../../../include/myNM.h"

int main(int argc, char* argv[])
{
	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif


	//+++++++++++++++++++
	//		Part1
	//+++++++++++++++++++
	
	printf("\n**************************************************");
	printf("\n|                     PART 1.                    |");
	printf("\n**************************************************\n");
	
	//Initialization
	double at[21] = { 0, 0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6, 3.8, 4 };
	double apos[21] = { -5.87,- 4.23,- 2.55,- 0.89,0.67,2.09,3.31,4.31,5.06,5.55,5.78,5.77,5.52,5.08,4.46,3.72,2.88,2.00,1.10,0.23, - 0.59 };
	double df[21] = { 0, };
	int m = 21;

	Matrix t = arr2Mat(at, 21, 1);
	Matrix pos = arr2Mat(apos, 21, 1);



	//primary code
	Matrix vel = gradient(t, pos);
	Matrix acc = gradient(t, vel);
	gradient1D(at, apos, df, m);
	Matrix mdf = arr2Mat(df, m, 1);



	//print result
	printMat(t, "t");
	printMat(pos, "pos");
	printMat(vel, "vel");
	printMat(acc, "acc");
	printMat(mdf, "result of gradient1D");






	//+++++++++++++++++++
	//		Part2
	//+++++++++++++++++++

	printf("\n**************************************************");
	printf("\n|                     PART 2.                    |");
	printf("\n**************************************************\n");

	//Initialization
	Matrix xin = arr2Mat(at, 21, 1);
	double epsilon = 1e-5;
	double x0 = 1;
	double x = 0;



	// primary code & print
	Matrix dydx = gradientFunc(myFunc, xin);
	printMat(xin, "xin");
	printMat(dydx, "dydx");

	printf("newtonRaphson result\n");
	x = newtonRaphsonFunc(myFunc, dmyFunc, x0, epsilon);


	system("pause");
	return 0;
}




