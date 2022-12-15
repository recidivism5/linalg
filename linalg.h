/*
Mat4 is column-major-order:
a0 b0 c0 d0 | 0 4 8  12
a1 b1 c1 d1 | 1 5 9  13
a2 b2 c2 d2 | 2 6 10 14
a3 b3 c3 d3 | 3 7 11 15

mat4Perspective and mat4LookAt assume a right-handed
coordinate system with camera looking down the
negative Z axis (OpenGL style):
	   y
	   |
	   |
	   |_ _ _ _ x
	   /
	  /
	 /
	z

*/
#define _USE_MATH_DEFINES
#include <math.h>

typedef union u_FVec3{
	struct {float x, y, z;};
	float arr[3];
}FVec3;
typedef union u_Quaternion{
	struct {float w, x, y, z;};
	float arr[4];
}Quaternion;
typedef union u_Mat4{
	struct {float 
		a0,a1,a2,a3,
		b0,b1,b2,b3,
		c0,c1,c2,c3,
		d0,d1,d2,d3;};
	float arr[16];
}Mat4;

float fvec3Dot(FVec3 a, FVec3 b);
float fvec3Length(FVec3 v);
FVec3 fvec3Cross(FVec3 a, FVec3 b);
FVec3 fvec3Add(FVec3 a, FVec3 b);
FVec3 fvec3Sub(FVec3 a, FVec3 b);
FVec3 fvec3Scale(FVec3 v, float s);
FVec3 fvec3Norm(FVec3 v);
FVec3 fvec3SetLength(FVec3 v, float length);
FVec3 fvec3Midpoint(FVec3 a, FVec3 b);
FVec3 fvec3Rotated(FVec3 v, Quaternion q);
float fvec3Dist(FVec3 a, FVec3 b);
float fvec3AngleBetween(FVec3 a, FVec3 b);
FVec3 clampEuler(FVec3 e);
Quaternion quatInverse(Quaternion q);
Quaternion quatNorm(Quaternion q);
Quaternion axisAngleToQuat(FVec3 axis, float angle);
Quaternion xToQuat(float angle);
Quaternion yToQuat(float angle);
Quaternion zToQuat(float angle);
Quaternion quatMul(Quaternion a, Quaternion b);
Quaternion eulerToQuat(FVec3 e);
Mat4 mat4Identity(void);
Mat4 mat4Basis(FVec3 i, FVec3 j, FVec3 k);
Mat4 quatToMat4(Quaternion q);
Mat4 mat4Transpose(Mat4 m);
Mat4 mat4Perspective(float fovRadians, float aspectRatio, float near, float far);
Mat4 mat4Scale(FVec3 v);
Mat4 mat4Pos(FVec3 v);
Mat4 xToMat4(float angle);
Mat4 yToMat4(float angle);
Mat4 zToMat4(float angle);
Mat4 mat4Mul(Mat4 a, Mat4 b);
Mat4 mat4LookAt(FVec3 eye, FVec3 target);