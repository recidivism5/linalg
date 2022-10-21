This library defines the following data types and functions:

```c
union FVec3 {
	struct { float x, y, z; };
	float arr[3];
};
struct Quaternion {
	float w, x, y, z;
};
void fvec3Normalize(union FVec3 *v);
void fvec3InvScale(union FVec3 s, union FVec3 *out);
void fvec3Scale(union FVec3 s, union FVec3 *out);
void fvec3Add(union FVec3 *dst, union FVec3 src);
union FVec3 fvec3_add(union FVec3 a, union FVec3 b);
union FVec3 fvec3_sub(union FVec3 a, union FVec3 b);
union FVec3 fvec3_scale(float s, union FVec3 v);
union FVec3 fvec3_norm(union FVec3 v);
union FVec3 fvec3_rescale(float s, union FVec3 v);
union FVec3 fvec3_midpoint(union FVec3 a, union FVec3 b);
float fvec3_length(union FVec3 v);
float fvec3_dist(union FVec3 a, union FVec3 b);
void clampEuler(union FVec3 *e);
void quatSetIdentity(struct Quaternion *q);
void quatInverse(struct Quaternion q, struct Quaternion *out);
void quatNormalize(struct Quaternion *q);
struct Quaternion quat_from(union FVec3 axis, float angle);
struct Quaternion quat_fromX(float angle);
struct Quaternion quat_fromY(float angle);
struct Quaternion quat_fromZ(float angle);
void quatMult(struct Quaternion a, struct Quaternion b, struct Quaternion *out);
void quatFromEuler(struct Quaternion *q, union FVec3 e);
struct Quaternion quat_from_euler(union FVec3 e);
void rotateFVec3(struct Quaternion q, union FVec3 *v);
union FVec3 fvec3_rotated(union FVec3 v, struct Quaternion q);
void quatToMat4(struct Quaternion q, float *m);
void mat4SetIdentity(float* m);
void mat4Transpose(float *m);
void mat4SetPerspective(float* m, float fov_radians, float aspect_ratio, float near, float far);
void mat4SetScale(float *m, float x, float y, float z);
void mat4SetPos(float* m, float x, float y, float z);
void mat4SetRotationX(float* m, float theta);
void mat4SetRotationY(float* m, float theta);
void mat4SetRotationZ(float* m, float theta);
void mat4Lookat(float *m, float x, float y, float z, float tx, float ty, float tz);
void mat4Mult(float* m0, float* m1, float* out);
```

- FVec3 is a vector of 3 floats
- Quaternion is a vector of 4 floats
- mat4 is an array of 16 floats, called "m" in functions using mat4s