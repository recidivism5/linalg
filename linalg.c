#include "linalg.h"

#define _USE_MATH_DEFINES
#include <math.h>

void fvec3Normalize(union FVec3 *v)
{
	float d = 1.0f / sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
	v->x *= d;
	v->y *= d;
	v->z *= d;
}
void fvec3InvScale(union FVec3 s, union FVec3 *out)
{
	out->x = 1.0f / s.x;
	out->y = 1.0f / s.y;
	out->z = 1.0f / s.z;
}
void fvec3Scale(union FVec3 s, union FVec3 *out)
{
	out->x *= s.x;
	out->y *= s.y;
	out->z *= s.z;
}
union FVec3 fvec3_add(union FVec3 a, union FVec3 b)
{
	union FVec3 r;
	r.x = a.x + b.x;
	r.y = a.y + b.y;
	r.z = a.z + b.z;
	return r;
}

void clampEuler(union FVec3 *e)
{
	if (e->x > 2*M_PI)
	{
		e->x -= 2*M_PI;
	}
	else if (e->x < -2*M_PI)
	{
		e->x += 2*M_PI;
	}

	if (e->y > 2*M_PI)
	{
		e->y -= 2*M_PI;
	}
	else if (e->y < -2*M_PI)
	{
		e->y += 2*M_PI;
	}

	if (e->z > 2*M_PI)
	{
		e->z -= 2*M_PI;
	}
	else if (e->z < -2*M_PI)
	{
		e->z += 2*M_PI;
	}
}

void quatSetIdentity(struct Quaternion *q)
{
	q->w = 1;
	q->x = 0;
	q->y = 0;
	q->z = 0;
}
void quatInverse(struct Quaternion q, struct Quaternion *out)
{
	out->w = q.w;
	out->x = -q.x;
	out->y = -q.y;
	out->z = -q.z;
}
void quatNormalize(struct Quaternion *q)
{
	float d = 1.0f / sqrt(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
	q->w *= d;
	q->x *= d;
	q->y *= d;
	q->z *= d;
}
void quatFrom(struct Quaternion *q, union FVec3 axis, float angle)
{
	float s = sin(angle * 0.5f);

	q->w = cos(angle * 0.5f);
	q->x = s * axis.x;
	q->y = s * axis.y;
	q->z = s * axis.z;
}
void quatFromX(struct Quaternion *q, float angle)
{
	union FVec3 a;
	a.x = 1;
	a.y = 0;
	a.z = 0;
	quatFrom(q, a, angle);
}
void quatFromY(struct Quaternion *q, float angle)
{
	union FVec3 a;
	a.x = 0;
	a.y = 1;
	a.z = 0;
	quatFrom(q, a, angle);
}
void quatFromZ(struct Quaternion *q, float angle)
{
	union FVec3 a;
	a.x = 0;
	a.y = 0;
	a.z = 1;
	quatFrom(q, a, angle);
}
void quatMult(struct Quaternion a, struct Quaternion b, struct Quaternion *out)
{
    out->w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
    out->x = a.x*b.w + a.w*b.x + a.y*b.z - a.z*b.y;
    out->y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
    out->z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
}
void quatFromEuler(struct Quaternion *q, union FVec3 e)
{
	float cx = cos(e.x * 0.5);
    float sx = sin(e.x * 0.5);
    float cy = cos(e.y * 0.5);
    float sy = sin(e.y * 0.5);
    float cz = cos(e.z * 0.5);
    float sz = sin(e.z * 0.5);

	q->w = cx*cy*cz + sx*sy*sz;
	q->x = sx*cy*cz - cx*sy*sz;
	q->y = cx*sy*cz + sx*cy*sz;
	q->z = cx*cy*sz - sx*sy*cz;
}
void quatRotateFVec3(struct Quaternion q, union FVec3 *v)
{
	float x,y,z;
	x = 2*(q.y*q.w*v->z - q.z*q.w*v->y + q.y*q.x*v->y + q.z*q.x*v->z) + q.w*q.w*v->x + q.x*q.x*v->x - q.y*q.y*v->x - q.z*q.z*v->x;
	y = 2*(q.x*q.y*v->x + q.z*q.y*v->z + q.w*q.z*v->x - q.x*q.w*v->z) + q.w*q.w*v->y - q.x*q.x*v->y + q.y*q.y*v->y - q.z*q.z*v->y;
	z = 2*(q.x*q.z*v->x + q.y*q.z*v->y - q.w*q.y*v->x + q.w*q.x*v->y) + q.w*q.w*v->z - q.x*q.x*v->z - q.y*q.y*v->z + q.z*q.z*v->z;

	v->x = x;
	v->y = y;
	v->z = z;
}
void quatToMat4(struct Quaternion q, float *m)
{
	m[0]  = 1 - 2*(q.y*q.y + q.z*q.z);
	m[1]  = 2*(q.x*q.y + q.z*q.w);
	m[2]  = 2*(q.x*q.z - q.y*q.w);
	m[3]  = 0;
 
	m[4]  = 2*(q.x*q.y - q.z*q.w);
	m[5]  = 1 - 2*(q.x*q.x + q.z*q.z);
	m[6]  = 2*(q.y*q.z + q.x*q.w);
	m[7]  = 0;
 
	m[8]  = 2*(q.x*q.z + q.y*q.w);
	m[9]  = 2*(q.y*q.z - q.x*q.w);
	m[10] = 1 - 2*(q.x*q.x + q.y*q.y);
	m[11] = 0;

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;
}

void mat4SetIdentity(float* m)
{
	m[0]  = 1.0f;
	m[1]  = 0.0f;
	m[2]  = 0.0f;
	m[3]  = 0.0f;

	m[4]  = 0.0f;
	m[5]  = 1.0f;
	m[6]  = 0.0f;
	m[7]  = 0.0f;

	m[8]  = 0.0f;
	m[9]  = 0.0f;
	m[10] = 1.0f;
	m[11] = 0.0f;

	m[12] = 0.0f;
	m[13] = 0.0f;
	m[14] = 0.0f;
	m[15] = 1.0f;
}

void mat4SetEqual(float* target, float* source)
{
	for (int i = 0; i < 16; i++)
	{
		target[i] = source[i];
	}
}

void mat4Transpose(float *m)
{
	/*
		0 4  8 12
		1 5  9 13
		2 6 10 14
		3 7 11 15
	*/
	float a1 = m[1],
		  a2 = m[2],
		  a3 = m[3],
		  a6 = m[6],
		  a7 = m[7],
		  a11 = m[11];
	
	m[1] = m[4];
	m[2] = m[8];
	m[3] = m[12];
	m[6] = m[9];
	m[7] = m[13];
	m[11] = m[14];

	m[4] = a1;
	m[8] = a2;
	m[12] = a3;
	m[9] = a6;
	m[13] = a7;
	m[14] = a11;
}

void mat4SetPerspective(float* m, float fov_radians, float aspect_ratio, float near, float far)
{
	float S = 1.0f / tan(fov_radians / 2.0f);
	float d = near - far;

	m[0] = S / aspect_ratio;
	m[5] = S;
	m[10] = (far + near) / d;

	m[11] = -1.0f;
	m[14] = (2.0f * far * near) / d;
	m[15] = 0.0f;

	m[1] = 0.0f;
	m[2] = 0.0f;
	m[3] = 0.0f;
	m[4] = 0.0f;
	
	m[6] = 0.0f;
	m[7] = 0.0f;
	m[8] = 0.0f;
	m[9] = 0.0f;

	m[12] = 0.0f;
	m[13] = 0.0f;
}

void mat4SetScale(float *m, float x, float y, float z){
	mat4SetIdentity(m);
	m[0] = x;
	m[5] = y;
	m[10] = z;
}

void mat4SetPos(float* m, float x, float y, float z){
	mat4SetIdentity(m);
	m[12] = x;
	m[13] = y;
	m[14] = z;
}

void mat4SetRotationX(float* m, float theta)
{
	mat4SetIdentity(m);
	float c = cos(theta);
	float s = sin(theta);
	m[0] = 1;
	m[1] = 0;
	m[2] = 0;
	m[4] = 0;
	m[5] = c;
	m[6] = s;
	m[8] = 0;
	m[9] = -s;
	m[10] = c;
}
void mat4SetRotationY(float* m, float theta)
{
	mat4SetIdentity(m);
	float c = cos(theta);
	float s = sin(theta);
	m[0] = c;
	m[1] = 0;
	m[2] = -s;
	m[4] = 0;
	m[5] = 1;
	m[6] = 0;
	m[8] = s;
	m[9] = 0;
	m[10] = c;
}
void mat4SetRotationZ(float* m, float theta)
{
	mat4SetIdentity(m);
	float c = cos(theta);
	float s = sin(theta);
	m[0] = c;
	m[1] = s;
	m[2] = 0;
	m[4] = -s;
	m[5] = c;
	m[6] = 0;
	m[8] = 0;
	m[9] = 0;
	m[10] = 1;
}
void mat4Lookat(float *m, float x, float y, float z, float tx, float ty, float tz)
{
	float kx = tx-x,
		  ky = ty-y,
		  kz = ty-z;
	float d = 1.0f / sqrt(kx*kx + ky*ky + kz*kz);
	kx *= d;
	ky *= d;
	kz *= d;
	
	float ix = -kz,
		  iy = 0.0f,
		  iz = kx;
	d = 1.0f / sqrt(ix*ix + iy*iy + iz*iz);
	ix *= d;
	iy *= d;
	iz *= d;

	float jx = iy*kz - iz*ky,
		  jy = iz*kx - ix*kz,
		  jz = ix*ky - iy*kx;
 
	m[0]  = ix;
	m[1]  = jx;
	m[2]  = -kx;
	m[3]  = 0.0f;
	 
	m[4]  = iy;
	m[5]  = jy;
	m[6]  = -ky;
	m[7]  = 0.0f;
 
	m[8]  = iz;
	m[9]  = jz;
	m[10] = -kz;
	m[11] = 0.0f;

	m[12] = -(ix*x + iy*y + iz*z);
	m[13] = -(jx*x + jy*y + jz*z);
	m[14] = kx*x + ky*y + kz*z;
	m[15] = 1.0f;
}

void mat4Mult(float* m0, float* m1, float* out)
{
	out[0] = m0[0] * m1[0] + m0[4] * m1[1] + m0[8] * m1[2] + m0[12] * m1[3];
	out[1] = m0[1] * m1[0] + m0[5] * m1[1] + m0[9] * m1[2] + m0[13] * m1[3];
	out[2] = m0[2] * m1[0] + m0[6] * m1[1] + m0[10] * m1[2] + m0[14] * m1[3];
	out[3] = m0[3] * m1[0] + m0[7] * m1[1] + m0[11] * m1[2] + m0[15] * m1[3];
	out[4] = m0[0] * m1[4] + m0[4] * m1[5] + m0[8] * m1[6] + m0[12] * m1[7];
	out[5] = m0[1] * m1[4] + m0[5] * m1[5] + m0[9] * m1[6] + m0[13] * m1[7];
	out[6] = m0[2] * m1[4] + m0[6] * m1[5] + m0[10] * m1[6] + m0[14] * m1[7];
	out[7] = m0[3] * m1[4] + m0[7] * m1[5] + m0[11] * m1[6] + m0[15] * m1[7];
	out[8] = m0[0] * m1[8] + m0[4] * m1[9] + m0[8] * m1[10] + m0[12] * m1[11];
	out[9] = m0[1] * m1[8] + m0[5] * m1[9] + m0[9] * m1[10] + m0[13] * m1[11];
	out[10] = m0[2] * m1[8] + m0[6] * m1[9] + m0[10] * m1[10] + m0[14] * m1[11];
	out[11] = m0[3] * m1[8] + m0[7] * m1[9] + m0[11] * m1[10] + m0[15] * m1[11];
	out[12] = m0[0] * m1[12] + m0[4] * m1[13] + m0[8] * m1[14] + m0[12] * m1[15];
	out[13] = m0[1] * m1[12] + m0[5] * m1[13] + m0[9] * m1[14] + m0[13] * m1[15];
	out[14] = m0[2] * m1[12] + m0[6] * m1[13] + m0[10] * m1[14] + m0[14] * m1[15];
	out[15] = m0[3] * m1[12] + m0[7] * m1[13] + m0[11] * m1[14] + m0[15] * m1[15];
}