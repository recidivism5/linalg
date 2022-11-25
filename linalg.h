typedef union u_FVec3 {
	struct { float x, y, z; };
	float arr[3];
} FVec3;
typedef union u_Quaternion {
	struct { float w, x, y, z; };
	float arr[4];
} Quaternion;
void fvec3Normalize(FVec3 *v);
void fvec3InvScale(FVec3 s, FVec3 *out);
void fvec3Scale(FVec3 s, FVec3 *out);
void fvec3Add(FVec3 *dst, FVec3 src);
void rotateFVec3(Quaternion q, FVec3 *v);
FVec3 fvec3_add(FVec3 a, FVec3 b);
FVec3 fvec3_sub(FVec3 a, FVec3 b);
FVec3 fvec3_scale(float s, FVec3 v);
FVec3 fvec3_norm(FVec3 v);
FVec3 fvec3_rescale(float s, FVec3 v);
FVec3 fvec3_midpoint(FVec3 a, FVec3 b);
FVec3 fvec3_rotated(FVec3 v, Quaternion q);
float fvec3_length(FVec3 v);
float fvec3_dist(FVec3 a, FVec3 b);
float fvec3_angle_between(FVec3 a, FVec3 b);
void clampEuler(FVec3 *e);
void quatSetIdentity(Quaternion *q);
void quatInverse(Quaternion q, Quaternion *out);
void quatNormalize(Quaternion *q);
Quaternion quat_from(FVec3 axis, float angle);
Quaternion quat_fromX(float angle);
Quaternion quat_fromY(float angle);
Quaternion quat_fromZ(float angle);
void quatMul(Quaternion a, Quaternion b, Quaternion *out);
void quatFromEuler(Quaternion *q, FVec3 e);
Quaternion quat_from_euler(FVec3 e);
void quatToMat4(Quaternion q, float *m);
void mat4SetIdentity(float* m);
void mat4Transpose(float *m);
void mat4SetPerspective(float* m, float fov_radians, float aspect_ratio, float near, float far);
void mat4SetScale(float *m, float x, float y, float z);
void mat4SetPos(float* m, float x, float y, float z);
void mat4SetRotationX(float* m, float theta);
void mat4SetRotationY(float* m, float theta);
void mat4SetRotationZ(float* m, float theta);
void mat4Lookat(float *m, float x, float y, float z, float tx, float ty, float tz);
void mat4Mul(float* m0, float* m1, float* out);