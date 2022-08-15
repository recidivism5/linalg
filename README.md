This library has mat4 and quaternion functions. Here's the header:

```
union FVec3 {
    struct { float x, y, z; };
    float arr[3];
};
void fvec3Normalize(union FVec3 *v);
void fvec3InvScale(union FVec3 s, union FVec3 *out);
void fvec3Scale(union FVec3 s, union FVec3 *out);
union FVec3 fvec3_add(union FVec3 a, union FVec3 b);
void clampEuler(union FVec3 *e);

struct IVec3 {
    int x, y, z;
};

struct Quaternion {
    float w, x, y, z;
};
void quatSetIdentity(struct Quaternion *q);
void quatInverse(struct Quaternion q, struct Quaternion *out);
void quatNormalize(struct Quaternion *q);
void quatFrom(struct Quaternion *q, union FVec3 axis, float angle);
void quatFromX(struct Quaternion *q, float angle);
void quatFromY(struct Quaternion *q, float angle);
void quatFromZ(struct Quaternion *q, float angle);
void quatMult(struct Quaternion a, struct Quaternion b, struct Quaternion *out);
void quatFromEuler(struct Quaternion *q, union FVec3 e);
void quatRotateFVec3(struct Quaternion q, union FVec3 *v);
void quatToMat4(struct Quaternion q, float *m);

void mat4SetIdentity(float* m);
void mat4SetEqual(float* target, float* source);
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