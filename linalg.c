#include "linalg.h"
#define FOR(var,count) for(int var = 0; var < (count); var++)

float fvec3Dot(FVec3 a, FVec3 b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
float fvec3Length(FVec3 v){
    return sqrtf(fvec3Dot(v,v));
}
FVec3 fvec3Cross(FVec3 a, FVec3 b){
    return (FVec3){a.y*b.z - a.z*b.y,
                   a.z*b.x - a.x*b.z,
                   a.x*b.y - a.y*b.x};
}
FVec3 fvec3Add(FVec3 a, FVec3 b){
    FOR(i,3) a.arr[i] += b.arr[i];
    return a;
}
FVec3 fvec3Sub(FVec3 a, FVec3 b){
    FOR(i,3) a.arr[i] -= b.arr[i];
    return a;
}
FVec3 fvec3Scale(FVec3 v, float s){
    FOR(i,3) v.arr[i] *= s;
    return v;
}
FVec3 fvec3Norm(FVec3 v){
    return fvec3Scale(v, 1.0f/fvec3Length(v));
}
FVec3 fvec3SetLength(FVec3 v, float length){
    return fvec3Scale(fvec3Norm(v), length);
}
FVec3 fvec3Midpoint(FVec3 a, FVec3 b){
    return fvec3Add(fvec3Scale(fvec3Sub(b,a),0.5f),a);
}
FVec3 fvec3Rotated(FVec3 v, Quaternion q){
    float ww = q.w*q.w,
          xx = q.x*q.x,
          yy = q.y*q.y,
          zz = q.z*q.z,

          wx = q.w*q.x,
          wy = q.w*q.y,
          wz = q.w*q.z,

          xy = q.x*q.y,
          xz = q.x*q.z,

          yz = q.y*q.z;
    FVec3 r;
    r.x = 2*(wy*v.z - wz*v.y + xy*v.y + xz*v.z) + ww*v.x + xx*v.x - yy*v.x - zz*v.x;
    r.y = 2*(xy*v.x + yz*v.z + wz*v.x - wx*v.z) + ww*v.y - xx*v.y + yy*v.y - zz*v.y;
    r.z = 2*(xz*v.x + yz*v.y - wy*v.x + wx*v.y) + ww*v.z - xx*v.z - yy*v.z + zz*v.z;
    return r;
}
float fvec3Dist(FVec3 a, FVec3 b){
    return fvec3Length(fvec3Sub(a,b));
}
float fvec3AngleBetween(FVec3 a, FVec3 b){
    return acosf(fvec3Dot(a,b)/(fvec3Length(a)*fvec3Length(b)));
}
FVec3 clampEuler(FVec3 e){
    float fp = 4*M_PI;
    FOR(i,3) if (e.arr[i] > fp) e.arr[i] -= fp;
           else if (e.arr[i] < -fp) e.arr[i] += fp;
    return e;
}
Quaternion quatInverse(Quaternion q){
    FOR(i,3) q.arr[1+i] = -q.arr[1+i];
    return q;
}
Quaternion quatNorm(Quaternion q){
    float d = 1.0f / sqrtf(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
    FOR(i,3) q.arr[i] *= d;
    return q;
}
Quaternion axisAngleToQuat(FVec3 axis, float angle){
    float s = sinf(angle*0.5f);
    Quaternion q;
    q.w = cosf(angle*0.5f);
    FOR(i,3) q.arr[1+i] = s*axis.arr[i];
    return q;
}
Quaternion xToQuat(float angle){
    return axisAngleToQuat((FVec3){1,0,0},angle);
}
Quaternion yToQuat(float angle){
    return axisAngleToQuat((FVec3){0,1,0},angle);
}
Quaternion zToQuat(float angle){
    return axisAngleToQuat((FVec3){0,0,1},angle);
}
Quaternion quatMul(Quaternion a, Quaternion b){
    Quaternion r;
    r.w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
    r.x = a.x*b.w + a.w*b.x + a.y*b.z - a.z*b.y;
    r.y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
    r.z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
    return r;
}
Quaternion eulerToQuat(FVec3 e){
    e = fvec3Scale(e,0.5f);
    float cx = cosf(e.x),
          sx = sinf(e.x),
          cy = cosf(e.y),
          sy = sinf(e.y),
          cz = cosf(e.z),
          sz = sinf(e.z);
    Quaternion q;
    q.w = cx*cy*cz + sx*sy*sz;
    q.x = sx*cy*cz - cx*sy*sz;
    q.y = cx*sy*cz + sx*cy*sz;
    q.z = cx*cy*sz - sx*sy*cz;
    return q;
}
Mat4 mat4Identity(void){
    return (Mat4){1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1};
}
Mat4 mat4Basis(FVec3 x, FVec3 y, FVec3 z){
    Mat4 m = mat4Identity();
    m.a0 = x.x; m.b0 = y.x; m.c0 = z.x;
    m.a1 = x.y; m.b1 = y.y; m.c1 = z.y;
    m.a2 = x.z; m.b2 = y.z; m.c2 = z.z;
    return m;
}
Mat4 quatToMat4(Quaternion q){
    Mat4 m = mat4Identity();
    m.a0 = 1 - 2*(q.y*q.y + q.z*q.z);
    m.a1 = 2*(q.x*q.y + q.z*q.w);
    m.a2 = 2*(q.x*q.z - q.y*q.w);

    m.b0 = 2*(q.x*q.y - q.z*q.w);
    m.b1 = 1 - 2*(q.x*q.x + q.z*q.z);
    m.b2 = 2*(q.y*q.z + q.x*q.w);

    m.c0 = 2*(q.x*q.z + q.y*q.w);
    m.c1 = 2*(q.y*q.z - q.x*q.w);
    m.c2 = 1 - 2*(q.x*q.x + q.y*q.y);
    return m;
}
Mat4 mat4Transpose(Mat4 m){
    float a1 = m.a1, a2 = m.a2, a3 = m.a3,
          b2 = m.b2, b3 = m.b3,
          c3 = m.c3;
    m.a1 = m.b0; m.a2 = m.c0; m.a3 = m.d0;
    m.b2 = m.c1; m.b3 = m.d1;
    m.c3 = m.d2;

    m.b0 = a1; m.c0 = a2; m.d0 = a3;
    m.c1 = b2; m.d1 = b3;
    m.d2 = c3;
    return m;
}
Mat4 mat4Perspective(float fovRadians, float aspectRatio, float near, float far){
    float s = 1.0f / tanf(fovRadians * 0.5f);
    float d = near - far;
    Mat4 m = mat4Identity();
    m.a0 = s / aspectRatio;
    m.b1 = s;
    m.c2 = (far + near) / d;
    m.d3 = 0;
    m.c3 = -1;
    m.d2 = (2.0f * far * near) / d;
    return m;
}
Mat4 mat4Scale(FVec3 v){
    Mat4 m = mat4Identity();
    m.a0 = v.x;
    m.b1 = v.y;
    m.c2 = v.z;
    return m;
}
Mat4 mat4Pos(FVec3 v){
    Mat4 m = mat4Identity();
    FOR(i,3) m.arr[3*4+i] = v.arr[i];
    return m;
}
Mat4 xToMat4(float angle){
    Mat4 m = mat4Identity();
    float s = sinf(angle), c = cosf(angle);
    m.b1 = c;
    m.b2 = s;
    m.c1 = -s;
    m.c2 = c;
    return m;
}
Mat4 yToMat4(float angle){
    Mat4 m = mat4Identity();
    float s = sinf(angle), c = cosf(angle);
    m.a0 = c;
    m.a2 = -s;
    m.c0 = s;
    m.c2 = c;
    return m;
}
Mat4 zToMat4(float angle){
    Mat4 m = mat4Identity();
    float s = sinf(angle), c = cosf(angle);
    m.a0 = c;
    m.a1 = s;
    m.b0 = -s;
    m.b1 = c;
    return m;
}
Mat4 mat4Mul(Mat4 a, Mat4 b){
    Mat4 m;
    FOR(i,4) FOR(j,4) m.arr[i*4+j] = a.arr[j]*b.arr[i*4] + a.arr[j+4]*b.arr[i*4+1] + a.arr[j+8]*b.arr[i*4+2] + a.arr[j+12]*b.arr[i*4+3];
    return m;
}
Mat4 mat4LookAt(FVec3 eye, FVec3 target){
    FVec3 z = fvec3Norm(fvec3Sub(eye,target)),
          x = fvec3Norm(fvec3Cross(z,(FVec3){0,1,0})),
          y = fvec3Cross(x,z);
    return mat4Mul(mat4Transpose(mat4Basis(x,y,z)),mat4Pos(fvec3Scale(eye,-1)));
}