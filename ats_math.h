#pragma once

#include <math.h>
#include <stdbool.h>

// ========================================== MACROS ========================================= //

#define PI (3.14159265359f)
#define TAU (6.28318530718f)

#define TO_RAD_MUL (0.01745329251f)
#define TO_DEG_MUL (57.2957795131f)

#define ToRad(deg) ((deg) * TO_RAD_MUL)
#define ToDeg(rad) ((rad) * TO_DEG_MUL)

#define ClampMin(n, min)    ((n) < (min)? (min) : (n))
#define ClampMax(n, max)    ((n) > (max)? (max) : (n))
#define Clamp(n, min, max)  ClampMax(ClampMin(n, min), max)
#define Clamp01(n)          Clamp(n, 0, 1)
#define Min(a, b)           ((a) < (b)? (a) : (b))
#define Max(a, b)           ((a) > (b)? (a) : (b))
#define Lerp(a, b, t)       ((a) + (float)(t) * ((b) - (a)))
#define Sign(n)             ((n) == 0? 0 : ((n) < 0? -1 : 1))

// =========================================== TYPES ========================================= //

typedef union {
  struct { float x, y; };
  float e[2];
} v2;

typedef union {
  struct { float x, y, z; };
  struct { float r, g, b; };
  struct { v2 xy; };
  float e[3];
} v3;

typedef union {
  struct { float x, y, z, w; };
  struct { float r, g, b, a; };
  struct { v3 rgb; };
  struct { v3 xyz; };
  float e[4];
} v4;

typedef union {
  struct { int x, y; };
  int e[2];
} v2i;

typedef union {
  struct { int x, y, z; };
  struct { v2i xy; };
  int e[3];
} v3i;

typedef union {
  struct { int x, y, z, w; };
  int e[4];
} v4i;

typedef union {
  struct { v2 x, y; };
  float e[4];
} m2;

typedef union {
  struct { v3 x, y, z; };
  float e[9];
} m3;

typedef union {
  struct { v4 x, y, z, w; };
  float e[16];
} m4;

typedef struct {
  v2 min;
  v2 max;
} r2;

typedef struct {
  v3 min;
  v3 max;
} r3;

typedef struct {
  v2i min;
  v2i max;
} r2i;

typedef struct {
  v3i min;
  v3i max;
} r3i;

typedef struct quat {
  float x, y, z, w;
} quat_t;

typedef struct circle {
  v2 pos;
  float rad;
} circle_t;

typedef struct sphere {
  v3 pos;
  float rad;
} sphere_t;

// ========================================= FUNCTIONS ====================================== //

#ifndef __cplusplus // C stuff

#define V2(...) ((v2) { __VA_ARGS__ })
#define V3(...) ((v3) { __VA_ARGS__ })
#define V4(...) ((v4) { __VA_ARGS__ })

#define V2i(...) ((v2i) { __VA_ARGS__ })
#define V3i(...) ((v3i) { __VA_ARGS__ })
#define V4i(...) ((v4i) { __VA_ARGS__ })

#define Quat(...) ((quat_t) { __VA_ARGS__ })

#define R2(...) ((r2) { __VA_ARGS__ })
#define R3(...) ((r3) { __VA_ARGS__ })

#define R2i(...) ((r2i) { __VA_ARGS__ })
#define R3i(...) ((r3i) { __VA_ARGS__ })

#define Circle(...) ((circle_t) { __VA_ARGS__ })
#define Sphere(...) ((sphere_t) { __VA_ARGS__ })

#define M2(...) ((m2) { __VA_ARGS__ })
#define M3(...) ((m3) { __VA_ARGS__ })
#define M4(...) ((m4) { __VA_ARGS__ })

#else // __cplusplus

static inline v2 V2(float x, float y)               { return { x, y }; }
static inline v3 V3(float x, float y, float z)        { return { x, y, z }; }
static inline v4 V4(float x, float y, float z, float w) { return { x, y, z, w }; }

static inline v2i V2i(int x, int y)               { return { x, y }; }
static inline v3i V3i(int x, int y, int z)        { return { x, y, z }; }
static inline v4i V4i(int x, int y, int z, int w) { return { x, y, z, w }; }

static inline quat_t Quat(float x, float y, float z, float w) { return { x, y, z, w }; }

static inline r2 R2(float ax, float ay, float bx, float by) { return { { ax, ay }, { bx, by } }; }
static inline r2 R2(v2 a, v2 b) { return { a, b }; }

static inline r3 R3(float ax, float ay, float az, float bx, float by, float bz) { return { { ax, ay, az }, { bx, by, bz } }; }
static inline r3 R3(v3 a, v3 b) { return { a, b }; }

static inline r2i R2i(int ax, int ay, int bx, int by) { return { { ax, ay }, { bx, by } }; }
static inline r2i R2i(v2i a, v2i b) { return { a, b }; }

static inline r3i R3i(int ax, int ay, int az, int bx, int by, int bz) { return { { ax, ay, az }, { bx, by, bz } }; }
static inline r3i R3i(v3i a, v3i b)  { return { a, b }; }

static inline circle_t Circle(float x, float y, float rad) { return { x, y, rad }; }
static inline circle_t Circle(v2 p, float rad) { return { p, rad }; }

static inline sphere_t Sphere(float x, float y, float z, float rad) { return { x, y, z, rad }; }
static inline sphere_t Sphere(v3 p, float rad) { return { p, rad }; }

#define M2(...) (m2 { __VA_ARGS__ })
#define M3(...) (m3 { __VA_ARGS__ })
#define M4(...) (m4 { __VA_ARGS__ })

#endif

static m2 m2_identity(void) {
  return M2(
    1, 0,
    0, 1);
}

static m3 m3_identity(void) {
  return M3(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1);
}

static m4 m4_identity(void) {
  return M4(
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1);
}

static quat_t quat_identity(void) {
  return Quat(0, 0, 0, 1);
}

static float sqrt32(float n) {
  float x = n * 0.5f;
  float y = n;
  int i = *(int*)&y;

  i = 0x5f3759df - (i >> 1);
  y = *(float*)&i;
  y = y * (1.5f - (x * y * y));

  return n * y;
}

static float rsqrt32(float n) {
  float x2 = n * 0.5f;
  float y  = n;
  int i  = *(long*)&y;           // evil floating point bit level hacking

  i = 0x5f3759df - (i >> 1);     // what the fuck? 
  y = *(float*) &i;
  y = y * (1.5f - (x2 * y * y)); // 1st iteration

  return y;
}

static float shortest_angle_distance(float a, float b) {
  float max = 2.0f * PI;
  float da  = fmodf(b - a, max);

  return fmodf(2.0f * da, max) - da;
}

static float lerp_angle(float a, float b, float t) {
  return a + shortest_angle_distance(a, b) * t;
}

static float sine_ease_in(float t) {
  return 1 - cosf((t * PI) / 2);
}

static float sine_ease_out(float t) {
  return sinf((t * PI) / 2);
}

static float sine_ease_in_out(float t) {
  return -0.5 * (cosf(PI * t) - 1);
}

static float quad_ease_in(float t) {
  return t * t;
}

static float quad_ease_out(float t) {
  return 1 - (1 - t) * (1 - t);
}

static float quad_ease_in_out(float t) {
  float k = -2 * t + 2;
  return (t < 0.5)? (2 * t * t) : (1 - 0.5 * k * k);
}

static float cubic_ease_in(float t) {
  return t * t * t;
}

static float cubic_ease_out(float t) {
  float k = 1 - t;
  return 1 - k * k * k;
}

static float cubic_ease_in_out(float t) {
  float k = -2 * t + 2;
  return (t < 0.5)? (4 * t * t * t) : (1 - 0.5 * k * k * k);
}

static float quart_ease_in(float t) {
  return t * t * t * t;
}

static float quart_ease_out(float t) {
  float k = 1 - t; 
  return 1 - k * k * k * k;
}

static float quart_ease_in_out(float t) {
  float k = -2 * t + 2;
  return (t < 0.5)? (8 * t * t * t * t) : (1 - 0.5 * k * k * k * k);
}

static float quint_ease_in(float t) {
  return t * t * t * t * t;
}

static float quint_ease_out(float t) {
  float k = 1 - t;
  return 1 - k * k * k * k * k;
}

static float quint_ease_in_out(float t) {
  float k = -2 * t + 2;
  return (t < 0.5)? (16 * t * t * t * t * t) : (1 - 0.5 * k * k * k * k * k);
}

static float expo_ease_in(float t) {
  return (t == 0)? 0 : powf(2, 10 * t - 10);
}

static float expo_ease_out(float t) {
  return (t == 1)? 1 : (1 - powf(2, -10 * t));
}

static float expo_ease_in_out(float t) {
  return (t == 0)? 0 : (t == 1)? 1 : t < 0.5? powf(2, 20 * t - 10) / 2 : (2 - powf(2, -20 * t + 10)) / 2;
}

static float circ_ease_in(float t) {
  return 1 - sqrt32(1 - (t * t));
}

static float circ_ease_out(float t) {
  return sqrt32(1 - (t - 1) * (t - 1));
}

static float circ_ease_in_out(float t) {
  float k = 2 * t;
  float l = -2 * t + 2;
  return (t < 0.5)? 0.5 * (1 - sqrt32(1 - k * k)) : 0.5 * (sqrt32(1 - l * l) + 1);
}

static float back_ease_in(float t) {
  float c1 = 1.70158;
  float c3 = c1 + 1;
  return c3 * t * t * t - c1 * t * t;
}

static float back_ease_out(float t) {
  float c1 = 1.70158;
  float c3 = c1 + 1;
  float k = t - 1;
  return 1 + c3 * k * k * k + c1 * k * k;
}

static float back_ease_in_out(float t) {
  float c1 = 1.70158;
  float c2 = c1 * 1.525;

  return (t < 0.5)?
    0.5 * (powf(2 * t, 2) * ((c2 + 1) * 2 * t - c2)) :
    0.5 * (pow(2 * t - 2, 2) * ((c2 + 1) * (t * 2 - 2) + c2) + 2);
}

static float elastic_ease_in(float t) {
  float c4 = (2 * PI) / 3;
  return (t == 0)?
    0 :
    (t == 1)?
    1 :
    -powf(2, 10 * t - 10) * sinf((t * 10 - 10.75) * c4);
}

static float elastic_ease_out(float t) {
  float c4 = (2 * PI) / 3;
  return t == 0?
    0 :
    t == 1?
    1 :
    powf(2, -10 * t) * sinf((t * 10 - 0.75) * c4) + 1;
}

static float elastic_ease_inout(float t) {
  float c5 = (2 * PI) / 4.5;
  return t == 0?
    0 :
    t == 1?
    1 :
    t < 0.5 ?
    -0.5 * (powf(2, 20 * t - 10)  * sinf((20 * t - 11.125) * c5)) :
    +0.5 * (powf(2, -20 * t + 10) * sinf((20 * t - 11.125) * c5)) + 1;
}

static float bounce_ease_out(float t) {
  float n1 = 7.5625;
  float d1 = 2.75;
  if (t < 1 / d1) {
    return n1 * t * t;
  } else if (t < 2 / d1) {
    t -= 1.5 / d1;
    return n1 * t * t + 0.75;
  } else if (t < 2.5 / d1) {
    t -= 2.25 / d1;
    return n1 * t * t + 0.9375;
  } else {
    t -= 2.625 / d1;
    return n1 * t * t + 0.984375;
  }
}

static float bounce_ease_in(float t) {
  return 1 - bounce_ease_out(t);
}

static float bounce_ease_in_out(float t) {
  return t < 0.5?
    0.5 * (1 - bounce_ease_out(1 - 2 * t)) :
    0.5 * (1 + bounce_ease_out(2 * t - 1));
}

// ---------- from array ---------- //

static v2 v2_from_array(const float* a) { return V2(a[0], a[1]); }
static v3 v3_from_array(const float* a) { return V3(a[0], a[1], a[2]); }
static v4 v4_from_array(const float* a) { return V4(a[0], a[1], a[2], a[3]); }

static v2i v2i_from_array(const int* a) { return V2i(a[0], a[1]); }
static v3i v3i_from_array(const int* a) { return V3i(a[0], a[1], a[2]); }
static v4i v4i_from_array(const int* a) { return V4i(a[0], a[1], a[2], a[3]); }

// ---------- unpack color ------------ //

static v3 v3_unpack_color(unsigned int color) {
  return V3(
    ((color & 0x000000ff) >> 0)  / 256.0f,
    ((color & 0x0000ff00) >> 8)  / 256.0f,
    ((color & 0x00ff0000) >> 16) / 256.0f);
}

static v4 v4_unpack_color(unsigned int color) {
  return V4(
    ((color & 0x000000ff) >> 0)  / 256.0f,
    ((color & 0x0000ff00) >> 8)  / 256.0f,
    ((color & 0x00ff0000) >> 16) / 256.0f,
    ((color & 0xff000000) >> 24) / 256.0f);
}

// --------- negate ---------- //

static v2 v2_neg(v2 u) { return V2(-u.x, -u.y); }
static v3 v3_neg(v3 u) { return V3(-u.x, -u.y, -u.z); }
static v4 v4_neg(v4 u) { return V4(-u.x, -u.y, -u.z, -u.w); }

static v2i v2i_neg(v2i u) { return V2i(-u.x, -u.y); }
static v3i v3i_neg(v3i u) { return V3i(-u.x, -u.y, -u.z); }
static v4i v4i_neg(v4i u) { return V4i(-u.x, -u.y, -u.z, -u.w); }

#ifdef __cplusplus

static v2 operator-(v2 u) { return { -u.x, -u.y }; }
static v3 operator-(v3 u) { return { -u.x, -u.y, -u.z }; }
static v4 operator-(v4 u) { return { -u.x, -u.y, -u.z, -u.w }; }

static v2i operator-(v2i u) { return { -u.x, -u.y }; }
static v3i operator-(v3i u) { return { -u.x, -u.y, -u.z }; }
static v4i operator-(v4i u) { return { -u.x, -u.y, -u.z, -u.w }; }

#endif

// ---------- addition ---------- //

static v2 v2_add(v2 a, v2 b) { return V2(a.x + b.x, a.y + b.y); }
static v3 v3_add(v3 a, v3 b) { return V3(a.x + b.x, a.y + b.y, a.z + b.z); }
static v4 v4_add(v4 a, v4 b) { return V4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w); }

static v2i v2i_add(v2i a, v2i b) { return V2i(a.x + b.x, a.y + b.y); }
static v3i v3i_add(v3i a, v3i b) { return V3i(a.x + b.x, a.y + b.y, a.z + b.z); }
static v4i v4i_add(v4i a, v4i b) { return V4i(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w); }

#ifdef __cplusplus

static v2 operator+(v2 a, v2 b) { return { a.x + b.x, a.y + b.y }; }
static v3 operator+(v3 a, v3 b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
static v4 operator+(v4 a, v4 b) { return { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w }; }

static v2i operator+(v2i a, v2i b) { return { a.x + b.x, a.y + b.y }; }
static v3i operator+(v3i a, v3i b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
static v4i operator+(v4i a, v4i b) { return { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w }; }

static v2& operator+=(v2& a, v2 b) { a = a + b; return a; }
static v3& operator+=(v3& a, v3 b) { a = a + b; return a; }
static v4& operator+=(v4& a, v4 b) { a = a + b; return a; }

static v2i& operator+=(v2i& a, v2i b) { a = a + b; return a; }
static v3i& operator+=(v3i& a, v3i b) { a = a + b; return a; }
static v4i& operator+=(v4i& a, v4i b) { a = a + b; return a; }

#endif

// -------- subtraction ------- //

static v2 v2_sub(v2 a, v2 b) { return V2(a.x - b.x, a.y - b.y); }
static v3 v3_sub(v3 a, v3 b) { return V3(a.x - b.x, a.y - b.y, a.z - b.z); }
static v4 v4_sub(v4 a, v4 b) { return V4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w); }

static v2i v2i_sub(v2i a, v2i b) { return V2i(a.x - b.x, a.y - b.y); }
static v3i v3i_sub(v3i a, v3i b) { return V3i(a.x - b.x, a.y - b.y, a.z - b.z); }
static v4i v4i_sub(v4i a, v4i b) { return V4i(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w); }

#ifdef __cplusplus

static v2 operator-(v2 a, v2 b) { return { a.x - b.x, a.y - b.y }; }
static v3 operator-(v3 a, v3 b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
static v4 operator-(v4 a, v4 b) { return { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w }; }

static v2i operator-(v2i a, v2i b) { return { a.x - b.x, a.y - b.y }; }
static v3i operator-(v3i a, v3i b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
static v4i operator-(v4i a, v4i b) { return { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w }; }

static v2& operator-=(v2& a, v2 b) { a = a - b; return a; }
static v3& operator-=(v3& a, v3 b) { a = a - b; return a; }
static v4& operator-=(v4& a, v4 b) { a = a - b; return a; }

static v2i& operator-=(v2i& a, v2i b) { a = a - b; return a; }
static v3i& operator-=(v3i& a, v3i b) { a = a - b; return a; }
static v4i& operator-=(v4i& a, v4i b) { a = a - b; return a; }

#endif

// -------- multiplication ------- //

static v2 v2_mul(v2 a, v2 b) { return V2(a.x * b.x, a.y * b.y); }
static v3 v3_mul(v3 a, v3 b) { return V3(a.x * b.x, a.y * b.y, a.z * b.z); }
static v4 v4_mul(v4 a, v4 b) { return V4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * a.w); }

static v2i v2i_mul(v2i a, v2i b) { return V2i(a.x * b.x, a.y * b.y); }
static v3i v3i_mul(v3i a, v3i b) { return V3i(a.x * b.x, a.y * b.y, a.z * a.z); }
static v4i v4i_mul(v4i a, v4i b) { return V4i(a.x * b.x, a.y * b.y, a.z * a.z, a.w * a.w); }

static v2 m2_mulv(m2 m, v2 u) {
  return V2(
    m.e[0] * u.x + m.e[2] * u.y,
    m.e[1] * u.x + m.e[3] * u.y);
}

static v3 m3_mulv(m3 m, v3 u) {
  return V3(
    m.e[0] * u.x + m.e[3] * u.y + m.e[6] * u.z,
    m.e[1] * u.x + m.e[4] * u.y + m.e[7] * u.z,
    m.e[2] * u.x + m.e[5] * u.y + m.e[8] * u.z);
}

static v4 m4_mulv(m4 m, v4 u) {
  return V4(
    m.e[0] * u.x + m.e[4] * u.y + m.e[8]  * u.z + m.e[12] * u.w,
    m.e[1] * u.x + m.e[5] * u.y + m.e[9]  * u.z + m.e[13] * u.w,
    m.e[2] * u.x + m.e[6] * u.y + m.e[10] * u.z + m.e[14] * u.w,
    m.e[3] * u.x + m.e[7] * u.y + m.e[11] * u.z + m.e[15] * u.w);
}

static m2 m2_mul(m2 a, m2 b) {
  return M2(
    a.e[0] * b.e[0] + a.e[2] * b.e[1],
    a.e[1] * b.e[0] + a.e[3] * b.e[1],
    a.e[0] * b.e[2] + a.e[2] * b.e[3],
    a.e[1] * b.e[2] + a.e[3] * b.e[3]);
}

static m3 m3_mul(m3 a, m3 b) {
  return M3(
    a.e[0] * b.e[0] + a.e[3] * b.e[1]  + a.e[6] * b.e[2],
    a.e[1] * b.e[0] + a.e[4] * b.e[1]  + a.e[7] * b.e[2],
    a.e[2] * b.e[0] + a.e[5] * b.e[1]  + a.e[8] * b.e[2],

    a.e[0] * b.e[3] + a.e[3] * b.e[4]  + a.e[6] * b.e[5],
    a.e[1] * b.e[3] + a.e[4] * b.e[4]  + a.e[7] * b.e[5],
    a.e[2] * b.e[3] + a.e[5] * b.e[4]  + a.e[8] * b.e[5],

    a.e[0] * b.e[6] + a.e[3] * b.e[7]  + a.e[6] * b.e[8],
    a.e[1] * b.e[6] + a.e[4] * b.e[7]  + a.e[7] * b.e[8],
    a.e[2] * b.e[6] + a.e[5] * b.e[7]  + a.e[8] * b.e[8]);
}

static m4 m4_mul(m4 a, m4 b) {
  return M4(
    a.e[0] * b.e[0]  + a.e[4] * b.e[1]  + a.e[8]  * b.e[2]  + a.e[12] * b.e[3],
    a.e[1] * b.e[0]  + a.e[5] * b.e[1]  + a.e[9]  * b.e[2]  + a.e[13] * b.e[3],
    a.e[2] * b.e[0]  + a.e[6] * b.e[1]  + a.e[10] * b.e[2]  + a.e[14] * b.e[3],
    a.e[3] * b.e[0]  + a.e[7] * b.e[1]  + a.e[11] * b.e[2]  + a.e[15] * b.e[3],

    a.e[0] * b.e[4]  + a.e[4] * b.e[5]  + a.e[8]  * b.e[6]  + a.e[12] * b.e[7],
    a.e[1] * b.e[4]  + a.e[5] * b.e[5]  + a.e[9]  * b.e[6]  + a.e[13] * b.e[7],
    a.e[2] * b.e[4]  + a.e[6] * b.e[5]  + a.e[10] * b.e[6]  + a.e[14] * b.e[7],
    a.e[3] * b.e[4]  + a.e[7] * b.e[5]  + a.e[11] * b.e[6]  + a.e[15] * b.e[7],

    a.e[0] * b.e[8]  + a.e[4] * b.e[9]  + a.e[8]  * b.e[10] + a.e[12] * b.e[11],
    a.e[1] * b.e[8]  + a.e[5] * b.e[9]  + a.e[9]  * b.e[10] + a.e[13] * b.e[11],
    a.e[2] * b.e[8]  + a.e[6] * b.e[9]  + a.e[10] * b.e[10] + a.e[14] * b.e[11],
    a.e[3] * b.e[8]  + a.e[7] * b.e[9]  + a.e[11] * b.e[10] + a.e[15] * b.e[11],

    a.e[0] * b.e[12] + a.e[4] * b.e[13] + a.e[8]  * b.e[14] + a.e[12] * b.e[15],
    a.e[1] * b.e[12] + a.e[5] * b.e[13] + a.e[9]  * b.e[14] + a.e[13] * b.e[15],
    a.e[2] * b.e[12] + a.e[6] * b.e[13] + a.e[10] * b.e[14] + a.e[14] * b.e[15],
    a.e[3] * b.e[12] + a.e[7] * b.e[13] + a.e[11] * b.e[14] + a.e[15] * b.e[15]);
}

static quat_t quat_mul(quat_t a, quat_t b) {
  return Quat(
    a.y * b.z - a.z * b.y + a.w * b.x + b.w * a.x,
    a.z * b.x - a.x * b.z + a.w * b.y + b.w * a.y,
    a.x * b.y - a.y * b.x + a.w * b.z + b.w * a.z,
    a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z);
}

#ifdef __cplusplus

static v2 operator*(v2 a, v2 b) { return { a.x * b.x, a.y * b.y }; }
static v3 operator*(v3 a, v3 b) { return { a.x * b.x, a.y * b.y, a.z * b.z }; }
static v4 operator*(v4 a, v4 b) { return { a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w }; }

static v2i operator*(v2i a, v2i b) { return { a.x * b.x, a.y * b.y }; }
static v3i operator*(v3i a, v3i b) { return { a.x * b.x, a.y * b.y, a.z * b.z }; }
static v4i operator*(v4i a, v4i b) { return { a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w }; }

static v2& operator*=(v2& a, v2 b) { a = a * b; return a; }
static v3& operator*=(v3& a, v3 b) { a = a * b; return a; }
static v4& operator*=(v4& a, v4 b) { a = a * b; return a; }

static v2i& operator*=(v2i& a, v2i b) { a = a * b; return a; }
static v3i& operator*=(v3i& a, v3i b) { a = a * b; return a; }
static v4i& operator*=(v4i& a, v4i b) { a = a * b; return a; }

static v2 operator*(m2 m, v2 u) { return m2_mulv(m, u); }
static v3 operator*(m3 m, v3 u) { return m3_mulv(m, u); }
static v4 operator*(m4 m, v4 u) { return m4_mulv(m, u); }

static m2 operator*(m2 a, m2 b) { return m2_mul(a, b); }
static m3 operator*(m3 a, m3 b) { return m3_mul(a, b); }
static m4 operator*(m4 a, m4 b) { return m4_mul(a, b); }

static quat_t operator*(quat_t a, quat_t b) { return quat_mul(a, b); }

#endif

// ------------ divition ------------ //

static v2 v2_div(v2 a, v2 b) { return V2(a.x / b.x, a.y / b.y); }
static v3 v3_div(v3 a, v3 b) { return V3(a.x / b.x, a.y / b.y, a.z / b.z); }
static v4 v4_div(v4 a, v4 b) { return V4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w); }

static v2i v2i_div(v2i a, v2i b) { return V2i(a.x / b.x, a.y / b.y); }
static v3i v3i_div(v3i a, v3i b) { return V3i(a.x / b.x, a.y / b.y, a.z / b.z); }
static v4i v4i_div(v4i a, v4i b) { return V4i(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w); }

#ifdef __cplusplus

static v2 operator/(v2 a, v2 b) { return { a.x / b.x, a.y / b.y }; }
static v3 operator/(v3 a, v3 b) { return { a.x / b.x, a.y / b.y, a.z / b.z }; }
static v4 operator/(v4 a, v4 b) { return { a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w }; }

static v2i operator/(v2i a, v2i b) { return { a.x / b.x, a.y / b.y }; }
static v3i operator/(v3i a, v3i b) { return { a.x / b.x, a.y / b.y, a.z / b.z }; }
static v4i operator/(v4i a, v4i b) { return { a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w }; }

static v2& operator/=(v2& a, v2 b) { a = a / b; return a; }
static v3& operator/=(v3& a, v3 b) { a = a / b; return a; }
static v4& operator/=(v4& a, v4 b) { a = a / b; return a; }

static v2i& operator/=(v2i& a, v2i b) { a = a / b; return a; }
static v3i& operator/=(v3i& a, v3i b) { a = a / b; return a; }
static v4i& operator/=(v4i& a, v4i b) { a = a / b; return a; }

#endif

// ------------- scaling ------------- //

static v2 v2_scale(v2 a, float s) { return V2(a.x * s, a.y * s); }
static v3 v3_scale(v3 a, float s) { return V3(a.x * s, a.y * s, a.z * s); }
static v4 v4_scale(v4 a, float s) { return V4(a.x * s, a.y * s, a.z * s, a.w * s); }

static v2i v2i_scale(v2i a, int s) { return V2i(a.x * s, a.y * s); }
static v3i v3i_scale(v3i a, int s) { return V3i(a.x * s, a.y * s, a.z * s); }
static v4i v4i_scale(v4i a, int s) { return V4i(a.x * s, a.y * s, a.z * s, a.w * s); }

#ifdef __cplusplus

static v2 operator*(v2 a, float s) { return { a.x * s, a.y * s }; }
static v3 operator*(v3 a, float s) { return { a.x * s, a.y * s, a.z * s }; }
static v4 operator*(v4 a, float s) { return { a.x * s, a.y * s, a.z * s, a.w * s }; }

static v2 operator*(float s, v2 a) { return { a.x * s, a.y * s }; }
static v3 operator*(float s, v3 a) { return { a.x * s, a.y * s, a.z * s }; }
static v4 operator*(float s, v4 a) { return { a.x * s, a.y * s, a.z * s, a.w * s }; }

static v2i operator*(v2i a, int s) { return { a.x * s, a.y * s }; }
static v3i operator*(v3i a, int s) { return { a.x * s, a.y * s, a.z * s }; }
static v4i operator*(v4i a, int s) { return { a.x * s, a.y * s, a.z * s, a.w * s }; }

static v2i operator*(int s, v2i a) { return { a.x * s, a.y * s }; }
static v3i operator*(int s, v3i a) { return { a.x * s, a.y * s, a.z * s }; }
static v4i operator*(int s, v4i a) { return { a.x * s, a.y * s, a.z * s, a.w * s }; }

static v2& operator*=(v2& a, float s) { a = a * s; return a; }
static v3& operator*=(v3& a, float s) { a = a * s; return a; }
static v4& operator*=(v4& a, float s) { a = a * s; return a; }

static v2i& operator*=(v2i& a, float s) { a = a * s; return a; }
static v3i& operator*=(v3i& a, float s) { a = a * s; return a; }
static v4i& operator*=(v4i& a, float s) { a = a * s; return a; }

#endif

// ----------- eq ------------ //

static bool v2i_eq(v2i a, v2i b) { return a.x == b.x && a.y == b.y; }
static bool v3i_eq(v3i a, v3i b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
static bool v4i_eq(v4i a, v4i b) { return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w; }

static bool v2i_neq(v2i a, v2i b) { return a.x != b.x || a.y != b.y; }
static bool v3i_neq(v3i a, v3i b) { return a.x != b.x || a.y != b.y || a.z != b.z; }
static bool v4i_neq(v4i a, v4i b) { return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w; }

#ifdef __cplusplus

static bool operator==(v2i a, v2i b) { return a.x == b.x && a.y == b.y; }
static bool operator==(v3i a, v3i b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
static bool operator==(v4i a, v4i b) { return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w; }

static bool operator!=(v2i a, v2i b) { return a.x != b.x || a.y != b.y; }
static bool operator!=(v3i a, v3i b) { return a.x != b.x || a.y != b.y || a.z != b.z; }
static bool operator!=(v4i a, v4i b) { return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w; }

#endif

// ----------- dot product ----------- //

static float v2_dot(v2 a, v2 b) { return a.x * b.x + a.y * b.y; }
static float v3_dot(v3 a, v3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static float v4_dot(v4 a, v4 b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

static int v2i_dot(v2i a, v2i b) { return a.x * b.x + a.y * b.y; }
static int v3i_dot(v3i a, v3i b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static int v4i_dot(v4i a, v4i b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

#ifdef __cplusplus

static float dot(v2 a, v2 b) { return a.x * b.x + a.y * b.y; }
static float dot(v3 a, v3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static float dot(v4 a, v4 b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

static int dot(v2i a, v2i b) { return a.x * b.x + a.y * b.y; }
static int dot(v3i a, v3i b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static int dot(v4i a, v4i b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

#endif

// ----------- length squared ----------- //

static float v2_len_sq(v2 u) { return v2_dot(u, u); }
static float v3_len_sq(v3 u) { return v3_dot(u, u); }
static float v4_len_sq(v4 u) { return v4_dot(u, u); }

static int v2i_len_sq(v2i u) { return v2i_dot(u, u); }
static int v3i_len_sq(v3i u) { return v3i_dot(u, u); }
static int v4i_len_sq(v4i u) { return v4i_dot(u, u); }

#ifdef __cplusplus

static float len_sq(v2 u) { return v2_dot(u, u); }
static float len_sq(v3 u) { return v3_dot(u, u); }
static float len_sq(v4 u) { return v4_dot(u, u); }

static int len_sq(v2i u) { return v2i_dot(u, u); }
static int len_sq(v3i u) { return v3i_dot(u, u); }
static int len_sq(v4i u) { return v4i_dot(u, u); }

#endif

// -------------- length -------------- //

static float v2_len(v2 u) { return sqrt32(v2_len_sq(u)); }
static float v3_len(v3 u) { return sqrt32(v3_len_sq(u)); }
static float v4_len(v4 u) { return sqrt32(v4_len_sq(u)); }

static float v2i_len(v2i u) { return sqrt32(v2i_len_sq(u)); }
static float v3i_len(v3i u) { return sqrt32(v3i_len_sq(u)); }
static float v4i_len(v4i u) { return sqrt32(v4i_len_sq(u)); }

#ifdef __cplusplus

static float len(v2 u) { return sqrt32(v2_len_sq(u)); }
static float len(v3 u) { return sqrt32(v3_len_sq(u)); }
static float len(v4 u) { return sqrt32(v4_len_sq(u)); }

static float len(v2i u) { return sqrt32(v2i_len_sq(u)); }
static float len(v3i u) { return sqrt32(v3i_len_sq(u)); }
static float len(v4i u) { return sqrt32(v4i_len_sq(u)); }

#endif

// -------------- distance squared -------------- //

static float v2_dist_sq(v2 a, v2 b) { return v2_len_sq(v2_sub(a, b)); }
static float v3_dist_sq(v3 a, v3 b) { return v3_len_sq(v3_sub(a, b)); }
static float v4_dist_sq(v4 a, v4 b) { return v4_len_sq(v4_sub(a, b)); }

static int v2i_dist_sq(v2i a, v2i b) { return v2i_len_sq(v2i_sub(a, b)); }
static int v3i_dist_sq(v3i a, v3i b) { return v3i_len_sq(v3i_sub(a, b)); }
static int v4i_dist_sq(v4i a, v4i b) { return v4i_len_sq(v4i_sub(a, b)); }

#ifdef __cplusplus

static float dist_sq(v2 a, v2 b) { return v2_len_sq(v2_sub(a, b)); }
static float dist_sq(v3 a, v3 b) { return v3_len_sq(v3_sub(a, b)); }
static float dist_sq(v4 a, v4 b) { return v4_len_sq(v4_sub(a, b)); }

static int dist_sq(v2i a, v2i b) { return v2i_len_sq(v2i_sub(a, b)); }
static int dist_sq(v3i a, v3i b) { return v3i_len_sq(v3i_sub(a, b)); }
static int dist_sq(v4i a, v4i b) { return v4i_len_sq(v4i_sub(a, b)); }

#endif

// ------------------ distance ------------------- //

static float v2_dist(v2 a, v2 b) { return sqrt32(v2_dist_sq(a, b)); }
static float v3_dist(v3 a, v3 b) { return sqrt32(v3_dist_sq(a, b)); }
static float v4_dist(v4 a, v4 b) { return sqrt32(v4_dist_sq(a, b)); }

#ifdef __cplusplus

static float dist(v2 a, v2 b) { return sqrt32(v2_dist_sq(a, b)); }
static float dist(v3 a, v3 b) { return sqrt32(v3_dist_sq(a, b)); }
static float dist(v4 a, v4 b) { return sqrt32(v4_dist_sq(a, b)); }

#endif

// -------------- manhattan distance -------------- //

static int v2i_manhattan(v2i a, v2i b) {
  v2i diff = v2i_sub(a, b);
  return (0x7ffffffff & diff.x) + (0x7ffffffff & diff.y);
}

static int v3i_manhattan(v3i a, v3i b) {
  v3i diff = v3i_sub(a, b);
  return (0x7ffffffff & diff.x) + (0x7ffffffff & diff.y) + (0x7ffffffff & diff.z);
}

#ifdef __cplusplus

static int manhattan(v2i a, v2i b) { return v2i_manhattan(a, b); }
static int manhattan(v3i a, v3i b) { return v3i_manhattan(a, b); }

#endif

// -------------- normalize --------------- //

static v2 v2_norm(v2 u) { return v2_scale(u, rsqrt32(v2_dot(u, u))); }
static v3 v3_norm(v3 u) { return v3_scale(u, rsqrt32(v3_dot(u, u))); }
static v4 v4_norm(v4 u) { return v4_scale(u, rsqrt32(v4_dot(u, u))); }

#ifdef __cplusplus

static v2 norm(v2 u) { return v2_scale(u, rsqrt32(v2_dot(u, u))); }
static v3 norm(v3 u) { return v3_scale(u, rsqrt32(v3_dot(u, u))); }
static v4 norm(v4 u) { return v4_scale(u, rsqrt32(v4_dot(u, u))); }

#endif

// -------------- project --------------- //

static v2 v2_project(v2 a, v2 b) {
  float d = v2_dot(b, b);
  if (d == 0) return V2(0);
  return v2_scale(b, v2_dot(a, b) / d);
}

static v3 v3_project(v3 a, v3 b) {
  float d = v3_dot(b, b);
  if (d == 0) return V3(0);
  return v3_scale(b, v3_dot(a, b) / d);
}

#ifdef __cplusplus

struct v2 project(v2 a, v2 b) { return v2_project(a, b); }
struct v3 project(v3 a, v3 b) { return v3_project(a, b); }

#endif

// ---------------- min ----------------- //

static v2
v2_min(v2 a, v2 b) {
  return V2(
    a.x < b.x? a.x : b.x,
    a.y < b.y? a.y : b.y);
}

static v3
v3_min(v3 a, v3 b) {
  return V3(
    a.x < b.x? a.x : b.x,
    a.y < b.y? a.y : b.y,
    a.z < b.z? a.z : b.z);
}

static v4
v4_min(v4 a, v4 b) {
  return V4(
    a.x < b.x? a.x : b.x,
    a.y < b.y? a.y : b.y,
    a.z < b.z? a.z : b.z,
    a.w < b.w? a.w : b.w);
}

static v2i
v2i_min(v2i a, v2i b) {
  return V2i(
    a.x < b.x? a.x : b.x,
    a.y < b.y? a.y : b.y);
}

static v3i
v3i_min(v3i a, v3i b) {
  return V3i(
    a.x < b.x? a.x : b.x,
    a.y < b.y? a.y : b.y,
    a.z < b.z? a.z : b.z);
}

static v4i
v4i_min(v4i a, v4i b) {
  return V4i(
    a.x < b.x? a.x : b.x,
    a.y < b.y? a.y : b.y,
    a.z < b.z? a.z : b.z,
    a.w < b.w? a.w : b.w);
}

#ifdef __cplusplus

static v2 min(v2 a, v2 b) { return v2_min(a, b); }
static v3 min(v3 a, v3 b) { return v3_min(a, b); }
static v4 min(v4 a, v4 b) { return v4_min(a, b); }

static v2i min(v2i a, v2i b) { return v2i_min(a, b); }
static v3i min(v3i a, v3i b) { return v3i_min(a, b); }
static v4i min(v4i a, v4i b) { return v4i_min(a, b); }

#endif

// ---------------- max ----------------- //

static v2
v2_max(v2 a, v2 b) {
  return V2(
    a.x > b.x? a.x : b.x,
    a.y > b.y? a.y : b.y);
}

static v3
v3_max(v3 a, v3 b) {
  return V3(
    a.x > b.x? a.x : b.x,
    a.y > b.y? a.y : b.y,
    a.z > b.z? a.z : b.z);
}

static v4
v4_max(v4 a, v4 b) {
  return V4(
    a.x > b.x? a.x : b.x,
    a.y > b.y? a.y : b.y,
    a.z > b.z? a.z : b.z,
    a.w > b.w? a.w : b.w);
}

static v2i
v2i_max(v2i a, v2i b) {
  return V2i(
    a.x > b.x? a.x : b.x,
    a.y > b.y? a.y : b.y);
}

static v3i
v3i_max(v3i a, v3i b) {
  return V3i(
    a.x > b.x? a.x : b.x,
    a.y > b.y? a.y : b.y,
    a.z > b.z? a.z : b.z);
}

static v4i
v4i_max(v4i a, v4i b) {
  return V4i(
    a.x > b.x? a.x : b.x,
    a.y > b.y? a.y : b.y,
    a.z > b.z? a.z : b.z,
    a.w > b.w? a.w : b.w);
}

#ifdef __cplusplus

static v2 max(v2 a, v2 b) { return v2_max(a, b); }
static v3 max(v3 a, v3 b) { return v3_max(a, b); }
static v4 max(v4 a, v4 b) { return v4_max(a, b); }

static v2i max(v2i a, v2i b) { return v2i_max(a, b); }
static v3i max(v3i a, v3i b) { return v3i_max(a, b); }
static v4i max(v4i a, v4i b) { return v4i_max(a, b); }

#endif

// ---------------- lerp ----------------- //

static v2 v2_lerp(v2 a, v2 b, float t) {
  return V2(
    a.x + t * (b.x - a.x),
    a.y + t * (b.y - a.y));
}

static v3 v3_lerp(v3 a, v3 b, float t) {
  return V3(
    a.x + t * (b.x - a.x),
    a.y + t * (b.y - a.y),
    a.z + t * (b.z - a.z));
}

static v4 v4_lerp(v4 a, v4 b, float t) {
  return V4(
    a.x + t * (b.x - a.x),
    a.y + t * (b.y - a.y),
    a.z + t * (b.z - a.z),
    a.w + t * (b.w - a.w));
}

#ifdef __cplusplus

static v2 lerp(v2 a, v2 b, float t) { return v2_lerp(a, b, t); }
static v3 lerp(v3 a, v3 b, float t) { return v3_lerp(a, b, t); }
static v4 lerp(v4 a, v4 b, float t) { return v4_lerp(a, b, t); }

#endif

// -------------- sign (-1, 0, 1) ------------------- //

static v2 v2_sign(v2 u) { return V2((float)Sign(u.x), (float)Sign(u.y)); }
static v3 v3_sign(v3 u) { return V3((float)Sign(u.x), (float)Sign(u.y), (float)Sign(u.z)); }
static v4 v4_sign(v4 u) { return V4((float)Sign(u.x), (float)Sign(u.y), (float)Sign(u.z), (float)Sign(u.w)); }

static v2i v2i_sign(v2i u) { return V2i(Sign(u.x), Sign(u.y)); }
static v3i v3i_sign(v3i u) { return V3i(Sign(u.x), Sign(u.y), Sign(u.z)); }
static v4i v4i_sign(v4i u) { return V4i(Sign(u.x), Sign(u.y), Sign(u.z), Sign(u.w)); }

#ifdef __cplusplus

static v2 sign(v2 u) { return v2_sign(u); }
static v3 sign(v3 u) { return v3_sign(u); }
static v4 sign(v4 u) { return v4_sign(u); }

static v2i sign(v2i u) { return v2i_sign(u); }
static v3i sign(v3i u) { return v3i_sign(u); }
static v4i sign(v4i u) { return v4i_sign(u); }

#endif

// --------------- cross ------------------- //

static v3 v3_cross(v3 a, v3 b) {
  return V3(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x);
}

#ifdef __cplusplus

static v3 cross(v3 a, v3 b) { return v3_cross(a, b); }

#endif

// --------------- get angle ------------------- //

static float v2_get_angle(v2 a, v2 b) {
  float det = a.x * b.y - b.x * a.y;
  float dot = a.x * b.x + a.y * b.y;
  return atan2f(det, dot);
}

#ifdef __cplusplus

static float get_angle(v2 a, v2 b) { return v2_get_angle(a, b); }

#endif

// ----------- keep min ---------- //

static v3 v3_keep_min(v3 u) {
  float dx = fabsf(u.x);
  float dy = fabsf(u.y);
  float dz = fabsf(u.z);
  if (dx <= dy && dx <= dz) return V3(u.x, 0, 0);
  if (dy <= dx && dy <= dz) return V3(0, u.y, 0);
  if (dz <= dx && dz <= dy) return V3(0, 0, u.z);
  return u;
}

#ifdef __cplusplus

static v3 keep_min(v3 u) { return v3_keep_min(u); }

#endif

// ----------- mask min ---------- //

static v3 v3_mask_min(v3 u) {
  float dx = fabsf(u.x);
  float dy = fabsf(u.y);
  float dz = fabsf(u.z);

  if (dx <= dy && dx <= dz) return V3(0, 1, 1);
  if (dy <= dx && dy <= dz) return V3(1, 0, 1);
  if (dz <= dx && dz <= dy) return V3(1, 1, 0);

  return V3(1, 1, 1);
}

#ifdef __cplusplus

static v3 mask_min(v3 u) { return v3_mask_min(u); }

#endif

// ------------------ transform/scale/rotate ------------------ //

static m2 m2_rotate(float angle) {
  float c = cosf(angle);
  float s = sinf(angle);

  return M2(c, s, -s, c);
}

static m3 m3_rotate(v3 axis, float angle) {
  float c = cosf(angle);
  float s = sinf(angle);
  float k = 1.0f - c;

  v3 sa = { s * axis.x, s * axis.y, s * axis.z };
  v3 omca = { k * axis.x, k * axis.y, k * axis.z };

  return M3(
    omca.x * axis.x + c,
    omca.x * axis.y - sa.z,
    omca.x * axis.z + sa.y,
    omca.y * axis.x + sa.z,
    omca.y * axis.y + c,
    omca.y * axis.z - sa.x,
    omca.z * axis.x - sa.y,
    omca.z * axis.y + sa.x,
    omca.z * axis.z + c);
}

static m4 m4_rotate(v3 axis, float angle) {
  float cosv = cosf(angle);
  float sinv = sinf(angle);
  float inv_cosv = 1.0f - cosv;

  v3 sa = { axis.x * sinv, axis.y * sinv, axis.z * sinv };
  v3 omca = { axis.x * inv_cosv, axis.y * inv_cosv, axis.z * inv_cosv };

  return M4(
    omca.x * axis.x + cosv,  omca.x * axis.y - sa.x,  omca.x * axis.z + sa.y, 0,
    omca.y * axis.x + sa.z,  omca.y * axis.y + cosv,  omca.y * axis.z - sa.x, 0,
    omca.z * axis.x - sa.y,  omca.z * axis.y + sa.x,  omca.z * axis.z + cosv, 0,
    0, 0, 0, 1);
}

static quat_t quat_rotate(v3 axis, float angle) {
  float s = sinf(0.5f * angle);
  v3 v = { s * axis.x, s * axis.y, s * axis.z };
  return Quat(v.x, v.y, v.z, cosf(0.5f * angle));
}

static m4 m4_translate(float x, float y, float z) {
  return M4(
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    x, y, z, 1);
}

static m4 m4_scale(float x, float y, float z) {
  return M4(
    x, 0, 0, 0,
    0, y, 0, 0,
    0, 0, z, 0,
    0, 0, 0, 1);
}

// --------------- from quat --------------- //

static m3 m3_from_quat(quat_t q) {
  float a = q.w;
  float b = q.x;
  float c = q.y;
  float d = q.z;

  float a2 = a * a;
  float b2 = b * b;
  float c2 = c * c;
  float d2 = d * d;

  return M3(
    a2 + b2 - c2 - d2,
    2.0f * (b * c + a * d),
    2.0f * (b * d - a * c),

    2.0f * (b * c - a * d),
    a2 - b2 + c2 - d2,
    2.0f * (c * d + a * b),

    2.0f * (b * d + a * c),
    2.0f * (c * d - a * b),
    a2 - b2 - c2 + d2);
}

static m4 m4_from_quat(quat_t q) {
  float a = q.w;
  float b = q.x;
  float c = q.y;
  float d = q.z;

  float a2 = a * a;
  float b2 = b * b;
  float c2 = c * c;
  float d2 = d * d;

  return M4(
    a2 + b2 - c2 - d2,
    2.0f * (b * c + a * d),
    2.0f * (b * d - a * c),
    0.0f,

    2.0f * (b * c - a * d),
    a2 - b2 + c2 - d2,
    2.0f * (c * d + a * b),
    0.0f,

    2.0f * (b * d + a * c),
    2.0f * (c * d - a * b),
    a2 - b2 - c2 + d2,
    0.0f,

    0.0f,
    0.0f,
    0.0f,
    1.0f);
}

// --------------- view matricies --------------- //

static m4 m4_ortho(float l, float r, float b, float t, float n, float f) {
  return M4(
    2 / (r - l), 0, 0, 0,
    0, 2 / (t - b), 0, 0,
    0, 0, -2 / (f - n), 0,
    -(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1);
}

static m4 m4_perspective(float y_fov, float aspect, float n, float f) {
  float a = 1.0f / tanf(y_fov / 2.0f);

  return M4(
    a / aspect, 0, 0, 0,
    0, a, 0, 0,
    0, 0, -((f + n) / (f - n)), -1,
    0, 0, -((2 * f * n) / (f - n)), 0);
}

static m4 m4_look_at(v3 eye, v3 center, v3 up) {
  v3 f = v3_norm(v3_sub(center, eye));
  v3 s = v3_norm(v3_cross(f, up));
  v3 t = v3_cross(s, f);
  m4 m;

  m.e[0]  =  s.x;
  m.e[1]  =  t.x;
  m.e[2]  = -f.x;
  m.e[3]  =  0;

  m.e[4]  =  s.y;
  m.e[5]  =  t.y;
  m.e[6]  = -f.y;
  m.e[7]  =  0;

  m.e[8]  =  s.z;
  m.e[9]  =  t.z;
  m.e[10] = -f.z;
  m.e[11] =  0;

  m.e[12] = -(m.e[0] * eye.x + m.e[4] * eye.y + m.e[8]  * eye.z);
  m.e[13] = -(m.e[1] * eye.x + m.e[5] * eye.y + m.e[9]  * eye.z);
  m.e[14] = -(m.e[2] * eye.x + m.e[6] * eye.y + m.e[10] * eye.z);
  m.e[15] = -(m.e[3] * eye.x + m.e[7] * eye.y + m.e[11] * eye.z - 1);

  return m;
}

// ----------------- plane/frustrum ------------------- //

typedef struct {
  float a;
  float b;
  float c;
  float d;
} plane_t;

typedef struct frustum {
  plane_t planes[6];
} frustum_t;

static plane_t plane_normalize(plane_t p) {
  float r_len = rsqrt32(p.a * p.a + p.b * p.b + p.c * p.c);

  p.a = p.a * r_len;
  p.b = p.b * r_len;
  p.c = p.c * r_len;
  p.d = p.d * r_len;

  return p;
}

static frustum_t frustum_create(m4 m) {
  frustum_t result;

  // left clipping plane
  result.planes[0].a = m.e[3]  + m.e[0];
  result.planes[0].b = m.e[7]  + m.e[4];
  result.planes[0].c = m.e[11] + m.e[8];
  result.planes[0].d = m.e[15] + m.e[12];

  // right clipping plane
  result.planes[1].a = m.e[3]  - m.e[0];
  result.planes[1].b = m.e[7]  - m.e[4];
  result.planes[1].c = m.e[11] - m.e[8];
  result.planes[1].d = m.e[15] - m.e[12];

  // top clipping plane
  result.planes[2].a = m.e[3]  - m.e[1];
  result.planes[2].b = m.e[7]  - m.e[5];
  result.planes[2].c = m.e[11] - m.e[9];
  result.planes[2].d = m.e[15] - m.e[13];

  // bottom clipping plane
  result.planes[3].a = m.e[3]  + m.e[1];
  result.planes[3].b = m.e[7]  + m.e[5];
  result.planes[3].c = m.e[11] + m.e[9];
  result.planes[3].d = m.e[15] + m.e[13];

  // near clipping plane
  result.planes[4].a = m.e[3]  + m.e[2];
  result.planes[4].b = m.e[7]  + m.e[6];
  result.planes[4].c = m.e[11] + m.e[10];
  result.planes[4].d = m.e[15] + m.e[14];

  // far clipping plane
  result.planes[5].a = m.e[3]  - m.e[2];
  result.planes[5].b = m.e[7]  - m.e[6];
  result.planes[5].c = m.e[11] - m.e[10];
  result.planes[5].d = m.e[15] - m.e[14];

  result.planes[0] = plane_normalize(result.planes[0]);
  result.planes[1] = plane_normalize(result.planes[1]);
  result.planes[2] = plane_normalize(result.planes[2]);
  result.planes[3] = plane_normalize(result.planes[3]);
  result.planes[4] = plane_normalize(result.planes[4]);
  result.planes[5] = plane_normalize(result.planes[5]);

  return result;
}

// ------------------ contains ------------------ //

static bool circle_contains(circle_t c, v2 pos) {
  float distance = v2_dist_sq(c.pos, pos);
  return distance < (c.rad * c.rad);
}

static bool sphere_contains(sphere_t s, v3 pos) {
  float distance = v3_dist_sq(s.pos, pos);
  return distance < (s.rad * s.rad);
}

static bool r2_contains(r2 rect, v2 pos) {
  if (pos.x < rect.min.x || pos.x > rect.max.x) return false;
  if (pos.y < rect.min.y || pos.y > rect.max.y) return false;
  return true;
}

static bool r3_contains(r3 rect, v3 pos) {
  if (pos.x < rect.min.x || pos.x > rect.max.x) return false;
  if (pos.y < rect.min.y || pos.y > rect.max.y) return false;
  if (pos.z < rect.min.z || pos.z > rect.max.z) return false;
  return true;
}

static bool r2i_contains(r2i rect, v2i pos) {
  if (pos.x < rect.min.x || pos.x > rect.max.x) return false;
  if (pos.y < rect.min.y || pos.y > rect.max.y) return false;
  return true;
}

static bool r3i_contains(r3i rect, v3i pos) {
  if (pos.x < rect.min.x || pos.x > rect.max.x) return false;
  if (pos.y < rect.min.y || pos.y > rect.max.y) return false;
  if (pos.z < rect.min.z || pos.z > rect.max.z) return false;
  return true;
}

static bool frustum_contains(frustum_t fs, v3 pos) {
  for (int i = 0; i < 6; i++) {
    if (fs.planes[i].a * pos.x + fs.planes[i].b * pos.y + fs.planes[i].c * pos.z + fs.planes[i].d <= 0)
      return false;
  }

  return true;
}

#ifdef __cplusplus

static bool contains(circle_t c, v2 pos) { return circle_contains(c, pos); }
static bool contains(sphere_t s, v3 pos) { return sphere_contains(s, pos); }
static bool contains(r2 rect, v2 pos)    { return r2_contains(rect, pos); }
static bool contains(r3 rect, v3 pos)    { return r3_contains(rect, pos); }
static bool contains(r2i rect, v2i pos)  { return r2i_contains(rect, pos); }
static bool contains(r3i rect, v3i pos)  { return r3i_contains(rect, pos); }
static bool contains(frustum_t fs, v3 pos) { return frustum_contains(fs, pos); }

#endif

// ------------------ intersect ------------------ //

static bool circle_intersect(circle_t a, circle_t b) {
  float dx  = b.pos.x - a.pos.x;
  float dy  = b.pos.y - a.pos.y;
  float rt  = a.rad + b.rad;
  return (dx * dx + dy * dy) < (rt * rt);
}

static bool sphere_intersect(sphere_t a, sphere_t b) {
  float dx = b.pos.x - a.pos.x;
  float dy = b.pos.y - a.pos.y;
  float dz = b.pos.z - a.pos.z;
  float rt = a.rad + b.rad;
  return (dx * dx + dy * dy) < (rt * rt);
}

static bool r2_intersect(r2 a, r2 b) {
  if (a.min.x > b.max.x || a.max.x < b.min.x) return false;
  if (a.min.y > b.max.y || a.max.y < b.min.y) return false;
  return true;
}

static bool r3_intersect(r3 a, r3 b) {
  if (a.min.x > b.max.x || a.max.x < b.min.x) return false;
  if (a.min.y > b.max.y || a.max.y < b.min.y) return false;
  if (a.min.z > b.max.z || a.max.z < b.min.z) return false;

  return true;
}

static bool r2i_intersect(r2i a, r2i b) {
  if (a.min.x > b.max.x || a.max.x < b.min.x) return false;
  if (a.min.y > b.max.y || a.max.y < b.min.y) return false;

  return true;
}

static bool r3i_intersect(r3i a, r3i b) {
  if (a.min.x > b.max.x || a.max.x < b.min.x) return false;
  if (a.min.y > b.max.y || a.max.y < b.min.y) return false;
  if (a.min.z > b.max.z || a.max.z < b.min.z) return false;

  return true;
}

static bool frustum_intersect_sphere(frustum_t fs, sphere_t sphere) {
  for (int i = 0; i < 6; i++) {
    if(fs.planes[i].a * sphere.pos.x + fs.planes[i].b * sphere.pos.y + fs.planes[i].c * sphere.pos.z + fs.planes[i].d <= -sphere.rad) {
      return false;
    }
  }
  return true;
}

static bool frustum_intersect_r3(frustum_t fs, r3 rect) {
  for (int i = 0; i < 6; i++) {
    if (fs.planes[i].a * rect.min.x + fs.planes[i].b * rect.min.y + fs.planes[i].c * rect.min.z + fs.planes[i].d > 0) continue;
    if (fs.planes[i].a * rect.max.x + fs.planes[i].b * rect.min.y + fs.planes[i].c * rect.min.z + fs.planes[i].d > 0) continue;
    if (fs.planes[i].a * rect.min.x + fs.planes[i].b * rect.max.y + fs.planes[i].c * rect.min.z + fs.planes[i].d > 0) continue;
    if (fs.planes[i].a * rect.max.x + fs.planes[i].b * rect.max.y + fs.planes[i].c * rect.min.z + fs.planes[i].d > 0) continue;
    if (fs.planes[i].a * rect.min.x + fs.planes[i].b * rect.min.y + fs.planes[i].c * rect.max.z + fs.planes[i].d > 0) continue;
    if (fs.planes[i].a * rect.max.x + fs.planes[i].b * rect.min.y + fs.planes[i].c * rect.max.z + fs.planes[i].d > 0) continue;
    if (fs.planes[i].a * rect.min.x + fs.planes[i].b * rect.max.y + fs.planes[i].c * rect.max.z + fs.planes[i].d > 0) continue;
    if (fs.planes[i].a * rect.max.x + fs.planes[i].b * rect.max.y + fs.planes[i].c * rect.max.z + fs.planes[i].d > 0) continue;
    return false;
  }
  return true;
}

#ifdef __cplusplus

static bool intersect(circle_t a, circle_t b)        { return circle_intersect(a, b); }
static bool intersect(sphere_t a, sphere_t b)        { return sphere_intersect(a, b); }
static bool intersect(r2 a, r2 b)                    { return r2_intersect(a, b); }
static bool intersect(r3 a, r3 b)                    { return r3_intersect(a, b); }
static bool intersect(r2i a, r2i b)                  { return r2i_intersect(a, b); }
static bool intersect(r3i a, r3i b)                  { return r3i_intersect(a, b); }
static bool intersect(frustum_t fs, sphere_t sphere) { return frustum_intersect_sphere(fs, sphere); }
static bool intersect(frustum_t fs, r3 rect)         { return frustum_intersect_r3(fs, rect); }

#endif

// ------------------- get overlap --------------- //

static r2 r2_get_overlap(r2 a, r2 b) {
  return R2(
    v2_max(a.min, b.min),
    v2_min(a.max, b.max));
}

static r3 r3_get_overlap(r3 a, r3 b) {
  return R3(
    v3_max(a.min, b.min),
    v3_min(a.max, b.max));
}

static r2i r2i_get_overlap(r2i a, r2i b) {
  return R2i(
    v2i_max(a.min, b.min),
    v2i_min(a.max, b.max));
}

static r3i r3i_get_overlap(r3i a, r3i b) {
  return R3i(
    v3i_max(a.min, b.min),
    v3i_min(a.max, b.max));
}

#ifdef __cplusplus

static r2  overlap(r2 a, r2 b)   { return r2_get_overlap(a, b); }
static r3  overlap(r3 a, r3 b)   { return r3_get_overlap(a, b); }
static r2i overlap(r2i a, r2i b) { return r2i_get_overlap(a, b); }
static r3i overlap(r3i a, r3i b) { return r3i_get_overlap(a, b); }

#endif

// -------------- get intersect vector ---------- //

static v2 circle_get_intersect_vector(circle_t a, circle_t b) {
  v2  delta = v2_sub(a.pos, b.pos);
  float depth = v2_len(delta) - (a.rad + b.rad);

  return v2_scale(delta, -depth);
}

static v3 sphere_get_intersect_vector(sphere_t a, sphere_t b) {
  v3 delta = v3_sub(a.pos, b.pos);
  float depth = v3_len(delta) - (a.rad + b.rad);

  return v3_scale(delta, -depth);
}

static v2 r2_get_intersect_vector(r2 a, r2 b) {
  r2 overlap  = r2_get_overlap(a, b);
  v2 delta = {
    0.5f * (a.min.x + a.max.x) - 0.5f * (b.min.x + b.max.x),
    0.5f * (a.min.y + a.max.y) - 0.5f * (b.min.y + b.max.y),
  };
  return v2_mul(v2_sub(overlap.max, overlap.min), v2_sign(delta));
}

static v3 r3_get_intersect_vector(r3 a, r3 b) {
  r3 overlap = r3_get_overlap(a, b);
  v3 delta = {
    0.5f * (a.min.x + a.max.x) - 0.5f * (b.min.x + b.max.x),
    0.5f * (a.min.y + a.max.y) - 0.5f * (b.min.y + b.max.y),
    0.5f * (a.min.z + a.max.z) - 0.5f * (b.min.z + b.max.z),
  };
  return v3_mul(v3_sub(overlap.max, overlap.min), v3_sign(delta));
}

static v2i r2i_get_intersect_vector(r2i a, r2i b) {
  r2i overlap = r2i_get_overlap(a, b);
  v2i delta = {
    (a.min.x + a.max.x) / 2 - (b.min.x + b.max.x) / 2,
    (a.min.y + a.max.y) / 2 - (b.min.y + b.max.y) / 2,
  };
  return v2i_mul(v2i_sub(overlap.max, overlap.min), v2i_sign(delta));
}

static v3i r3i_get_intersect_vector(r3i a, r3i b) {
  r3i overlap = r3i_get_overlap(a, b);
  v3i delta = {
    (a.min.x + a.max.x) / 2 - (b.min.x + b.max.x) / 2,
    (a.min.y + a.max.y) / 2 - (b.min.y + b.max.y) / 2,
    (a.min.z + a.max.z) / 2 - (b.min.z + b.max.z) / 2,
  };
  return v3i_mul(v3i_sub(overlap.max, overlap.min), v3i_sign(delta));
}

#ifdef __cplusplus

static v2  get_intersect_vector(circle_t a, circle_t b) { return circle_get_intersect_vector(a, b); }
static v3  get_intersect_vector(sphere_t a, sphere_t b) { return sphere_get_intersect_vector(a, b); }
static v2  get_intersect_vector(r2 a, r2 b)             { return r2_get_intersect_vector(a, b); }
static v3  get_intersect_vector(r3 a, r3 b)             { return r3_get_intersect_vector(a, b); }
static v2i get_intersect_vector(r2i a, r2i b)           { return r2i_get_intersect_vector(a, b); }
static v3i get_intersect_vector(r3i a, r3i b)           { return r3i_get_intersect_vector(a, b); }

#endif

// ---------------------- random ------------------------ //

static unsigned int rand_u32(unsigned int* state) {
  unsigned int x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  return *state = x;
}

// [min, max)
static int rand_i32(unsigned int* state, int min, int max) {
  return min + rand_u32(state) % (max - min);
}

static float rand_f32(unsigned int* state, float min, float max) {
  return min + ((float)rand_u32(state) / (float)0xffffffff) * (max - min); 
}

static v2 rand_unit_v2(unsigned int* state) {
  v2 out = { rand_f32(state, -1, 1), rand_f32(state, -1, 1) };
  return v2_norm(out);
}

static v3 rand_unit_v3(unsigned int* state) {
  v3 out = { rand_f32(state, -1, 1), rand_f32(state, -1, 1), rand_f32(state, -1, 1) };
  return v3_norm(out);
}

static v2 rand_v2(unsigned int* state, float min, float max) {
  return v2_scale(rand_unit_v2(state), rand_f32(state, min, max));
}

static v3 rand_v3(unsigned int* state, float min, float max) {
  return v3_scale(rand_unit_v3(state), rand_f32(state, min, max));
}

// ----------------------- hash ------------------------- //

static unsigned int hash_u32(unsigned int a) {
  a = (a ^ 61) ^ (a >> 16);
  a = a + (a << 3);
  a = a ^ (a >> 4);
  a = a * 0x27d4eb2d;
  a = a ^ (a >> 15);

  return a;
}

static unsigned int hash_i32(int a) {
  union { unsigned int u; int i; } convert;
  convert.i = a;
  return hash_u32(convert.u);
}

static unsigned int hash_str(const char* str) {
  unsigned int hash = 0;
  while (*str) {
    hash = (hash << 7) + (hash >> 25) + *str++;
  }
  return hash + (hash >> 16);
}

static const unsigned int crc_table[] = {
  0x00000000, 0x04c11db7, 0x09823b6e, 0x0d4326d9, 0x130476dc, 0x17c56b6b, 0x1a864db2, 0x1e475005,
  0x2608edb8, 0x22c9f00f, 0x2f8ad6d6, 0x2b4bcb61, 0x350c9b64, 0x31cd86d3, 0x3c8ea00a, 0x384fbdbd,
  0x4c11db70, 0x48d0c6c7, 0x4593e01e, 0x4152fda9, 0x5f15adac, 0x5bd4b01b, 0x569796c2, 0x52568b75,
  0x6a1936c8, 0x6ed82b7f, 0x639b0da6, 0x675a1011, 0x791d4014, 0x7ddc5da3, 0x709f7b7a, 0x745e66cd,
  0x9823b6e0, 0x9ce2ab57, 0x91a18d8e, 0x95609039, 0x8b27c03c, 0x8fe6dd8b, 0x82a5fb52, 0x8664e6e5,
  0xbe2b5b58, 0xbaea46ef, 0xb7a96036, 0xb3687d81, 0xad2f2d84, 0xa9ee3033, 0xa4ad16ea, 0xa06c0b5d,
  0xd4326d90, 0xd0f37027, 0xddb056fe, 0xd9714b49, 0xc7361b4c, 0xc3f706fb, 0xceb42022, 0xca753d95,
  0xf23a8028, 0xf6fb9d9f, 0xfbb8bb46, 0xff79a6f1, 0xe13ef6f4, 0xe5ffeb43, 0xe8bccd9a, 0xec7dd02d,
  0x34867077, 0x30476dc0, 0x3d044b19, 0x39c556ae, 0x278206ab, 0x23431b1c, 0x2e003dc5, 0x2ac12072,
  0x128e9dcf, 0x164f8078, 0x1b0ca6a1, 0x1fcdbb16, 0x018aeb13, 0x054bf6a4, 0x0808d07d, 0x0cc9cdca,
  0x7897ab07, 0x7c56b6b0, 0x71159069, 0x75d48dde, 0x6b93dddb, 0x6f52c06c, 0x6211e6b5, 0x66d0fb02,
  0x5e9f46bf, 0x5a5e5b08, 0x571d7dd1, 0x53dc6066, 0x4d9b3063, 0x495a2dd4, 0x44190b0d, 0x40d816ba,
  0xaca5c697, 0xa864db20, 0xa527fdf9, 0xa1e6e04e, 0xbfa1b04b, 0xbb60adfc, 0xb6238b25, 0xb2e29692,
  0x8aad2b2f, 0x8e6c3698, 0x832f1041, 0x87ee0df6, 0x99a95df3, 0x9d684044, 0x902b669d, 0x94ea7b2a,
  0xe0b41de7, 0xe4750050, 0xe9362689, 0xedf73b3e, 0xf3b06b3b, 0xf771768c, 0xfa325055, 0xfef34de2,
  0xc6bcf05f, 0xc27dede8, 0xcf3ecb31, 0xcbffd686, 0xd5b88683, 0xd1799b34, 0xdc3abded, 0xd8fba05a,
  0x690ce0ee, 0x6dcdfd59, 0x608edb80, 0x644fc637, 0x7a089632, 0x7ec98b85, 0x738aad5c, 0x774bb0eb,
  0x4f040d56, 0x4bc510e1, 0x46863638, 0x42472b8f, 0x5c007b8a, 0x58c1663d, 0x558240e4, 0x51435d53,
  0x251d3b9e, 0x21dc2629, 0x2c9f00f0, 0x285e1d47, 0x36194d42, 0x32d850f5, 0x3f9b762c, 0x3b5a6b9b,
  0x0315d626, 0x07d4cb91, 0x0a97ed48, 0x0e56f0ff, 0x1011a0fa, 0x14d0bd4d, 0x19939b94, 0x1d528623,
  0xf12f560e, 0xf5ee4bb9, 0xf8ad6d60, 0xfc6c70d7, 0xe22b20d2, 0xe6ea3d65, 0xeba91bbc, 0xef68060b,
  0xd727bbb6, 0xd3e6a601, 0xdea580d8, 0xda649d6f, 0xc423cd6a, 0xc0e2d0dd, 0xcda1f604, 0xc960ebb3,
  0xbd3e8d7e, 0xb9ff90c9, 0xb4bcb610, 0xb07daba7, 0xae3afba2, 0xaafbe615, 0xa7b8c0cc, 0xa379dd7b,
  0x9b3660c6, 0x9ff77d71, 0x92b45ba8, 0x9675461f, 0x8832161a, 0x8cf30bad, 0x81b02d74, 0x857130c3,
  0x5d8a9099, 0x594b8d2e, 0x5408abf7, 0x50c9b640, 0x4e8ee645, 0x4a4ffbf2, 0x470cdd2b, 0x43cdc09c,
  0x7b827d21, 0x7f436096, 0x7200464f, 0x76c15bf8, 0x68860bfd, 0x6c47164a, 0x61043093, 0x65c52d24,
  0x119b4be9, 0x155a565e, 0x18197087, 0x1cd86d30, 0x029f3d35, 0x065e2082, 0x0b1d065b, 0x0fdc1bec,
  0x3793a651, 0x3352bbe6, 0x3e119d3f, 0x3ad08088, 0x2497d08d, 0x2056cd3a, 0x2d15ebe3, 0x29d4f654,
  0xc5a92679, 0xc1683bce, 0xcc2b1d17, 0xc8ea00a0, 0xd6ad50a5, 0xd26c4d12, 0xdf2f6bcb, 0xdbee767c,
  0xe3a1cbc1, 0xe760d676, 0xea23f0af, 0xeee2ed18, 0xf0a5bd1d, 0xf464a0aa, 0xf9278673, 0xfde69bc4,
  0x89b8fd09, 0x8d79e0be, 0x803ac667, 0x84fbdbd0, 0x9abc8bd5, 0x9e7d9662, 0x933eb0bb, 0x97ffad0c,
  0xafb010b1, 0xab710d06, 0xa6322bdf, 0xa2f33668, 0xbcb4666d, 0xb8757bda, 0xb5365d03, 0xb1f740b4
};

static unsigned int hash_mem(const void *data, unsigned int size) {
  const unsigned char *d = (const unsigned char*)data;
  unsigned int crc = 0xFFFFFFFF;
  while (size--) {
    unsigned int index = (crc ^ *(d++)) & 0xFF;
    crc = (crc >> 8) ^ crc_table[index];
  }
  return crc ^ 0xFFFFFFFF;
}

#define HASH_PRIME0 3323784421u
#define HASH_PRIME1 1449091801u
#define HASH_PRIME2 4280703257u
#define HASH_PRIME3 1609059329u

static unsigned int hash_v2i(v2i k) {
  unsigned int a = hash_i32(k.x);
  unsigned int b = hash_i32(k.y);
  return (a * HASH_PRIME0) ^ (b * HASH_PRIME1);
}

static unsigned int hash_v3i(v3i k) {
  unsigned int a = hash_i32(k.x);
  unsigned int b = hash_i32(k.y);
  unsigned int c = hash_i32(k.z);
  return (a * HASH_PRIME0) ^ (b * HASH_PRIME1) ^ (c * HASH_PRIME2);
}

static unsigned int hash_v4i(v4i k) {
  unsigned int a = hash_i32(k.x);
  unsigned int b = hash_i32(k.y);
  unsigned int c = hash_i32(k.z);
  unsigned int d = hash_i32(k.w);
  return (a * HASH_PRIME0) ^ (b * HASH_PRIME1) ^ (c * HASH_PRIME2) ^ (d * HASH_PRIME3);
}

