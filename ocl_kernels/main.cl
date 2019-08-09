// FinalityLabs - 2019
// fixed 768 size prime-field arithmetic library (add, sub, mul, pow)
// Montgomery reduction parameters:
// B = 2^32 (Because our digits are uint32)


typedef uint uint32;
typedef ulong uint64;

typedef uint64 limb;
typedef uint64 limb2;

typedef struct {
  limb v[12];
} int768;

#define FIELD_LIMBS (12)
#define LIMB_BITS (64)
#define LIMB_MAX (0xfffffffffffffffful)

// Montgomery form of 1 = (1 * R mod P)
//
#define mnt4753_ONE ((int768){{0x98a8ecabd9dc6f42,0x91cd31c65a034686,0x97c3e4a0cd14572e,0x79589819c788b601,0xed269c942108976f,0x1e0f4d8acf031d68,0x320c3bb713338559,0x598b4302d2f00a62,0x4074c9cbfd8ca621,0xfa47edb3865e88c,0x95455fb31ff9a195,0x7b479ec8e242}})
//
//#define mnt6753_ONE ((int768){0x,0x,0x,0x,0x,0x,0x,0x,0x,0x,0x,0x})

#define mnt4753_ZERO ((int768){{0}})
#define mnt6753_ZERO ((int768){{0}})

#define mnt6753_INV_Fq ((ulong)0xc90776e23fffffff)
//#define mnt4753_INV_Fq ((uint){{0xe45e7fff}})
#define mnt4753_INV_Fq ((ulong)0xf2044cfbe45e7fff)
//#define mnt6753_INV_Fq ((uint){{0x3fffffff}})

#define mnt4753_Q ((int768){{0x5e9063de245e8001,0xe39d54522cdd119f,0x638810719ac425f0,0x685acce9767254a4,0xb80f0da5cb537e38,0xb117e776f218059d,0x99d124d9a15af79d,0x7fdb925e8a0ed8d,0x5eb7e8f96c97d873,0xb7f997505b8fafed,0x10229022eee2cdad,0x1c4c62d92c411}})
//#define mnt6753_Q ((int768){{0x40000001,0xd90776e2,0xfa13a4f,0x4ea09917,0x3f005797,0xd6c381bc,0x34993aa4,0xb9dff976,0x29212636,0x3eebca94,0xc859a99b,0xb26c5c28,0xa15af79d,0x99d124d9,0xe8a0ed8d,0x7fdb925,0x6c97d873,0x5eb7e8f9,0x5b8fafed,0xb7f99750,0xeee2cdad,0x10229022,0x2d92c411,0x1c4c6}})

limb bitreverse(limb n, limb bits) {
  limb r = 0;
  for(int i = 0; i < bits; i++) {
    r = (r << 1) | (n & 1);
    n >>= 1;
  }
  return r;
}


bool EXPONENT_get_bit(int768 l, uint i) {
  return (l.v[11 - i / 64] >> (63 - (i % 64))) & 1;
}

uint EXPONENT_get_bits(int768 l, uint skip, uint window) {
  uint ret = 0;
  for(uint i = 0; i < window; i++) {
    ret <<= 1;
    ret |= EXPONENT_get_bit(l, skip + i);
  }
  return ret;
}


// Greater than or equal
bool int768_gte(int768 a, int768 b) {
  for(int i = FIELD_LIMBS - 1; i >= 0; i--){
    if(a.v[i] > b.v[i])
      return true;
    if(a.v[i] < b.v[i])
      return false;
  }
  return true;
}

limb add_with_carry(limb a, limb *b) {
  limb lo, hi;
  asm("add.cc.u64 %0, %2, %3;\r\n"
      "addc.u64 %1, 0, 0;\r\n"
      : "=l"(lo), "=l"(hi) : "l"(a), "l"(*b));
  *b = hi;
  return lo;
}

limb add2_with_carry(limb a, limb b, bool *c) {
  limb lo, hi, cc = *c;
  asm("add.cc.u64 %0, %2, %3;\r\n"
      "addc.u64 %1, 0, 0;\r\n"
      "add.cc.u64 %0, %0, %4;\r\n"
      "addc.u64 %1, %1, 0;\r\n"
      : "=l"(lo), "=l"(hi) : "l"(a), "l"(b), "l"(cc));
  *c = hi;
  return lo;
}

limb mac_with_carry(limb a, limb b, limb c, limb *d) {
  limb lo, hi;
  asm("mad.lo.cc.u64 %0, %2, %3, %4;\r\n"
      "madc.hi.u64 %1, %2, %3, 0;\r\n"
      "add.cc.u64 %0, %0, %5;\r\n"
      "addc.u64 %1, %1, 0;\r\n"
      : "=l"(lo), "=l"(hi) : "l"(a), "l"(b), "l"(c), "l"(*d));
  *d = hi;
  return lo;
}


// Equals
bool int768_eq(int768 a, int768 b) {
  for(int i = 0; i < FIELD_LIMBS; i++)
    if(a.v[i] != b.v[i])
      return false;
  return true;
}

// Normal addition
int768 int768_add_(int768 a, int768 b) {
  bool carry = 0;
  for(int i = 0; i < FIELD_LIMBS; i++) {
    //limb2 sum = (limb2)a.v[i] + b.v[i] + carry;
    //a.v[i] = sum & LIMB_MAX;
    //carry = sum >> LIMB_BITS;
    //"this implementation fails for some inputs"
    limb old = a.v[i];
    a.v[i] += b.v[i] + carry;
    carry = carry ? old >= a.v[i] : old > a.v[i];
  }
  return a;
}

// Normal subtraction
int768 int768_sub_(int768 a, int768 b) {
  bool borrow = 0;
  for(int i = 0; i < FIELD_LIMBS; i++) {
    //limb2 sub = (limb2)a.v[i] - b.v[i] - borrow;
    //a.v[i] = sub & LIMB_MAX;
    //borrow = (sub >> LIMB_BITS) & 1;
    // "still works for sub but removing for consistency"
    limb old = a.v[i];
    a.v[i] -= b.v[i] + borrow;
    borrow = borrow ? old <= a.v[i] : old < a.v[i];
  }
  return a;
}

int768 int768_reduce4(ulong *limbs) {
  // Montgomery reduction
  bool carry2 = 0;
  for(uchar i = 0; i < FIELD_LIMBS; i++) {
    limb u = mnt4753_INV_Fq * limbs[i];
    limb carry = 0;
    for(uchar j = 0; j < FIELD_LIMBS; j++) {
      limbs[i + j] = mac_with_carry(u, mnt4753_Q.v[j], limbs[i + j], &carry);
    }
    limbs[i + FIELD_LIMBS] = add2_with_carry(limbs[i + FIELD_LIMBS], carry, &carry2);
  }

  // Divide by R
  int768 result;
  // this breaks amd compiler
  for(uchar i = 0; i < FIELD_LIMBS; i++) result.v[i] = limbs[i+FIELD_LIMBS];

  if(int768_gte(result, mnt4753_Q))
    result = int768_sub_(result, mnt4753_Q);

  return result;
}

//int768 int768_reduce6(ulong *limbs) {
//  // Montgomery reduction
//  bool carry2 = 0;
//  for(uchar i = 0; i < FIELD_LIMBS; i++) {
//    limb u = mnt6753_INV_Fq * limbs[i];
//    limb carry = 0;
//    for(uchar j = 0; j < FIELD_LIMBS; j++) {
//      limbs[i + j] = mac_with_carry(u, mnt6753_Q.v[j], limbs[i + j], &carry);
//    }
//    limbs[i + FIELD_LIMBS] = add2_with_carry(limbs[i + FIELD_LIMBS], carry, &carry2);
//  }

//  // Divide by R
//  int768 result;
//  // this breaks amd compiler
//  for(uchar i = 0; i < FIELD_LIMBS; i++) result.v[i] = limbs[i+FIELD_LIMBS];

//  if(int768_gte(result, mnt6753_Q))
//    result = int768_sub_(result, mnt6753_Q);

// return result;
//}

// Modular multiplication
int768 int768_mul4(int768 a, int768 b) {
  // Long multiplication
  limb res[FIELD_LIMBS * 2] = {0};
  for(uchar i = 0; i < FIELD_LIMBS; i++) {
    limb carry = 0;
    for(uchar j = 0; j < FIELD_LIMBS; j++)
      res[i + j] = mac_with_carry(a.v[i], b.v[j], res[i + j], &carry);
    res[i + FIELD_LIMBS] = carry;
  }

  return int768_reduce4(res);
}

// // Modular multiplication
//int768 int768_mul6(int768 a, int768 b) {
//  // Long multiplication
//  limb res[FIELD_LIMBS * 2] = {0};
//  for(uchar i = 0; i < FIELD_LIMBS; i++) {
//    limb carry = 0;
//    for(uchar j = 0; j < FIELD_LIMBS; j++)
//      res[i + j] = mac_with_carry(a.v[i], b.v[j], res[i + j], &carry);
//    res[i + FIELD_LIMBS] = carry;
// }
//
//  return int768_reduce6(res);
//}

// Modular negation
int768 int768_neg4(int768 a) {
  return int768_sub_(mnt4753_Q, a);
}

//int768 int768_neg6(int768 a) {
//  return int768_sub_(mnt6753_Q, a);
//}

// Modular subtraction
int768 int768_sub4(int768 a, int768 b) {
  int768 res = int768_sub_(a, b);
  if(!int768_gte(a, b)) res = int768_add_(res, mnt4753_Q);
  return res;
}

//int768 int768_sub6(int768 a, int768 b) {
//  int768 res = int768_sub_(a, b);
//  if(!int768_gte(a, b)) res = int768_add_(res, mnt6753_Q);
//  return res;
//}

// Modular addition
int768 int768_add4(int768 a, int768 b) {
  //return int768_sub4(a, int768_neg4(b));
  int768 tmp = int768_neg4(b);
  int768 res = int768_sub_(a, tmp);
  if(!int768_gte(a, tmp)) res = int768_add_(res, mnt4753_Q);
  return res;
}

//int768 int768_add6(int768 a, int768 b) {
//  //return int768_sub4(a, int768_neg6(b));
//  int768 tmp = int768_neg6(b);
//  int768 res = int768_sub_(a, tmp);
//  if(!int768_gte(a, tmp)) res = int768_add_(res, mnt6753_Q);
//  return res;
//}

// Modular exponentiation
int768 int768_pow(int768 base, uint32 exponent) {
  int768 res = mnt4753_ONE;
  while(exponent > 0) {
    if (exponent & 1)
      res = int768_mul4(res, base);
    exponent = exponent >> 1;
    base = int768_mul4(base, base);
  }
  return res;
}

int768 int768_pow_cached(__global int768 *bases, uint32 exponent) {
  int768 res = mnt4753_ONE;
  uint32 i = 0;
  while(exponent > 0) {
    if (exponent & 1)
      res = int768_mul4(res, bases[i]);
    exponent = exponent >> 1;
    i++;
  }
  return res;
}


// Fq2 MNT4753 arithmetics
//
typedef struct {
  int768 c0;
  int768 c1;
} Fq2;

#define Fq2_ZERO ((Fq2){mnt4753_ZERO, mnt4753_ZERO})
#define Fq2_ONE ((Fq2){mnt4753_ONE, FIELD_ZERO})

// Montgomery non_residue
#define non_residue ((int768){{0xa4e2d91fa3162657,0xbc938a1c0b935ff7,0x8a5a6ad599bbfb8a,0xf06f5292be9a4027,0xe2c8ca944b7535ff,0x737f39a7ace06d7a,0xbd2b99bf158cdead,0x74193bb2fc4dbe53,0x29c6846f9a5ce658,0xa36dab30ca7dbf57,0x641e2bafd304cb88,0xf450877b312e}});
// Integer non_residue
//#define non_residue ((int768){{0xD}});

bool Fq2_eq(Fq2 a, Fq2 b) {
  return 
  int768_eq(a.c0, b.c0) && int768_eq(a.c1, b.c1);
}

Fq2 Fq2_neg(Fq2 a) {
  a.c0 = int768_neg4(a.c0);
  a.c1 = int768_neg4(a.c1);
  return a;
}

Fq2 Fq2_sub(Fq2 a, Fq2 b) {
  a.c0 = int768_sub4(a.c0, b.c0);
  a.c1 = int768_sub4(a.c1, b.c1);
  return a;
}

Fq2 Fq2_add(Fq2 a, Fq2 b) {
  a.c0 = int768_add4(a.c0, b.c0);
  a.c1 = int768_add4(a.c1, b.c1);
  return a;
}

Fq2 Fq2_mul(Fq2 _a, Fq2 _b) {
  int768 residue = non_residue;
  int768 A = _b.c0;
  int768 a = _a.c0;
  int768 B = _b.c1;
  int768 b = _a.c1;

  int768 aA = int768_mul4(a, A);
  int768 bB = int768_mul4(b, B);

  Fq2 res = Fq2_ZERO;

  res.c0 = int768_add4(int768_mul4(_a.c0, _b.c0), int768_mul4(residue, int768_mul4(_a.c1, _b.c1)));
  
  res.c1 = int768_mul4(int768_add4(_a.c0, _a.c1), int768_add4(_b.c0, _b.c1));
  res.c1 = int768_sub4(res.c1, aA);
  res.c1 = int768_sub4(res.c1, bB);

  return res;
}


// Fq3 arithmetics
//

//typedef struct {
//  int768 c0;
//  int768 c1;
//  int768 c2;
//} Fq3;

//#define Fq3_ZERO ((Fq3){mnt6753_ZERO, mnt6753_ZERO, mnt6753_ZERO})
//#define Fq3_ONE ((Fq3){mnt6753_ONE, mnt6753_ZERO, mnt6753_ZERO})

// Montgomery non_residue
//#define non_residue ((int768){{0x4768931cfff9c7d4,0xc45e46d6ada96ca0,0x479b0bdb0b3c0107,0x362a089610f8d41b,0xdbafcec2c8a91aaf,0x78428b0ff9d96a06,0xf2e4472a9080c353,0xc9006ed33f0e971c,0x794d9d10bdb7288,0x3c1e44cab5419e2c,0x49b5fc6c81f4560c,0x1c287777c30ba}});
// Integer non_residue
//#define non_residue ((int768){{0xD}});


//Fq3 int768_fq3_add(Fq3 x, Fq3 y) {
//  Fq3 res = Fq3_ZERO;
//  res.c0 = int768_add6(x.c0, y.c0); 
//  res.c1 = int768_add6(x.c1, y.c1);
//  res.c2 = int768_add6(x.c2, y.c2);
//  return res;
//}

//Fq3 int768_fq3_mul(Fq3 a, Fq3 b) {
//  Fq3 res = Fq3_ZERO;
//  int768 residue = non_residue;
//
//  int768 a0_b0 = int768_mul6(a.c0, b.c0);
//  int768 a0_b1 = int768_mul6(a.c0, b.c1); 
//  int768 a0_b2 = int768_mul6(a.c0, b.c2);
//
//  int768 a1_b0 = int768_mul6(a.c1, b.c0);
//  int768 a1_b1 = int768_mul6(a.c1, b.c1); 
//  int768 a1_b2 = int768_mul6(a.c1, b.c2);
//
//  int768 a2_b0 = int768_mul6(a.c2, b.c0);
//  int768 a2_b1 = int768_mul6(a.c2, b.c1); 
//  int768 a2_b2 = int768_mul6(a.c2, b.c2);
//
//  res.c0 = int768_add6(a0_b0, int768_mul6(residue, int768_add6(a1_b2, a2_b1)));
//  res.c1 = int768_add6(a0_b1, int768_add6(a1_b0, int768_mul6(residue, a2_b2)));
//  res.c2 = int768_add6(a0_b2, int768_add6(a1_b1, a2_b0));
//  return res;
//}

// eliptic curve arithmetics
//

typedef struct {
  int768 X_;
  int768 Y_;
  int768 Z_;
} MNT_G1;

// affine coord zero
#define G1_ZERO ((MNT_G1){mnt4753_ZERO, mnt4753_ONE, mnt4753_ZERO})
#define G1_COEFF_A ((int768){{0x3151d957b3b8de84,0x239a638cb4068d0d,0x2f87c9419a28ae5d,0xf2b130338f116c03,0xda4d392842112ede,0x3c1e9b159e063ad1,0x6418776e26670ab2,0xb3168605a5e014c4,0x80e99397fb194c42,0x1f48fdb670cbd118,0x2a8abf663ff3432a,0xf68f3d91c485}})

bool is_zero(MNT_G1 a) {
  if(!int768_eq(a.X_, mnt4753_ZERO)) return false;
  // if(int768_eq(a.Y_, mnt4753_ONE)) return false;
  if(!int768_eq(a.Z_, mnt4753_ZERO)) return false; 
  return true;
}

// dont think we need this
bool G1_eq(MNT_G1 a, MNT_G1 b) {
  if(!int768_eq(a.X_, b.X_)) return false;
  if(!int768_eq(a.X_, b.X_)) return false;
  if(!int768_eq(a.X_, b.X_)) return false;
  return true;
}

MNT_G1 G1_add4(MNT_G1 a, MNT_G1 b) {
  if(is_zero(a)) return b;
  if(is_zero(b)) return a;

  MNT_G1 res = G1_ZERO;

  int768 X1_Z2 = int768_mul4(a.X_, b.Z_);
  int768 X2_Z1 = int768_mul4(a.Z_, b.X_);

  int768 Y1_Z2 = int768_mul4(a.Y_, b.Z_);
  int768 Y2_Z1 = int768_mul4(a.Z_, b.Y_);

  // double case 
  if(int768_eq(X1_Z2, X2_Z1) && int768_eq(Y1_Z2, Y2_Z1)) {
    int768 XX = int768_mul4(a.X_, a.X_); // todo special case squaring
    int768 ZZ = int768_mul4(a.Z_, a.Z_);
    int768 TXX = int768_add4(XX, XX);
    TXX = int768_add4(TXX, XX);
    int768 wz = int768_mul4(G1_COEFF_A, ZZ);
    int768 w = int768_add4(wz, TXX);
    int768 Y1_Z1 = int768_mul4(a.Y_, a.Z_);
    int768 s = int768_add4(Y1_Z1, Y1_Z1);
    int768 ss = int768_mul4(s, s);
    int768 sss = int768_mul4(s, ss);
    int768 R = int768_mul4(a.Y_, s);
    int768 RR = int768_mul4(R, R);
    int768 XR = int768_add4(a.X_, R);
    int768 XRXR = int768_mul4(XR, XR);
    XRXR = int768_sub4(XRXR, XX);
    int768 B = int768_sub4(XRXR, RR);
    int768 ww = int768_mul4(w, w);
    int768 BB = int768_add4(B, B);
    int768 h = int768_sub4(ww, BB);
    int768 X3 = int768_mul4(h, s);
    int768 b_h = int768_sub4(B, h);
    int768 wbh = int768_mul4(w, b_h);
    int768 RRRR = int768_add4(RR, RR);
    int768 Y3 = int768_sub4(wbh, RRRR);
    res.X_ = X3;
    res.Y_ = Y3;
    res.Z_ = sss;
    return res;
  }

  // add case
  int768 Z1_Z2 = int768_mul4(a.Z_, b.Z_);
  int768 u = int768_sub4(Y2_Z1, Y1_Z2);
  int768 uu = int768_mul4(u, u);
  int768 v = int768_sub4(X2_Z1, X1_Z2);
  int768 vv = int768_mul4(v,v);
  int768 vvv = int768_mul4(v,vv);
  int768 R = int768_mul4(vv, X1_Z2);
  int768 vvvR = int768_add4(vvv, R);
  vvvR = int768_add4(vvvR, R);
  int768 A = int768_sub4(int768_mul4(uu, Z1_Z2), vvvR);
  int768 X3 = int768_mul4(v, A);
  int768 vvvY1Z2 = int768_mul4(vvv, Y1_Z2);
  int768 Y3 = int768_sub4(int768_mul4(u, int768_sub4(R, A)), vvvY1Z2); 
  int768 Z3 = int768_mul4(vvv, Z1_Z2);

  res.X_ = X3;
  res.Y_ = Y3;
  res.Z_ = Z3;
  return res;
}

MNT_G1 G1_double4(MNT_G1 a) {
  if(int768_eq(a.Z_, mnt4753_ZERO)) return a;

  MNT_G1 res = G1_ZERO;

  int768 XX = int768_mul4(a.X_, a.X_); // todo special case squaring
  int768 ZZ = int768_mul4(a.Z_, a.Z_);
  int768 TXX = int768_add4(XX, XX);
  TXX = int768_add4(TXX, XX);
  int768 wz = int768_mul4(G1_COEFF_A, ZZ);
  int768 w = int768_add4(wz, TXX);
  int768 Y1_Z1 = int768_mul4(a.Y_, a.Z_);
  int768 s = int768_add4(Y1_Z1, Y1_Z1);
  int768 ss = int768_mul4(s, s);
  int768 sss = int768_mul4(s, ss);
  int768 R = int768_mul4(a.Y_, s);
  int768 RR = int768_mul4(R, R);
  int768 XR = int768_add4(a.X_, R);
  int768 XRXR = int768_mul4(XR, XR);
  XRXR = int768_sub4(XRXR, XX);
  int768 B = int768_sub4(XRXR, RR);
  int768 ww = int768_mul4(w, w);
  int768 BB = int768_add4(B, B);
  int768 h = int768_sub4(ww, BB);
  int768 X3 = int768_mul4(h, s);
  int768 b_h = int768_sub4(B, h);
  int768 wbh = int768_mul4(w, b_h);
  int768 RRRR = int768_add4(RR, RR);
  int768 Y3 = int768_sub4(wbh, RRRR);
  res.X_ = X3;
  res.Y_ = Y3;
  res.Z_ = sss;
  return res;
}


MNT_G1 G1_mixed_add4(MNT_G1 a, MNT_G1 b) {
  if(int768_eq(a.Z_, mnt4753_ZERO)) {
    a.X_ = b.X_;
    a.Y_ = b.Y_;
    a.Z_ = mnt4753_ONE;
    return a;
  }

  MNT_G1 res = G1_ZERO;
  int768 X1_Z2 = a.X_;
  int768 X2_Z1 = int768_mul4(a.Z_, b.X_);

  int768 Y1_Z2 = a.Y_;
  int768 Y2_Z1 = int768_mul4(a.Z_, b.Y_);

  if(int768_eq(X1_Z2, X2_Z1) && int768_eq(Y1_Z2, Y2_Z1)) {
    return G1_double4(a);
  }

  int768 u = int768_sub4(Y2_Z1, a.Y_);
  int768 uu = int768_mul4(u, u);
  int768 v = int768_sub4(X2_Z1, a.X_);
  int768 vv = int768_mul4(v,v);
  int768 vvv = int768_mul4(v,vv);
  int768 R = int768_mul4(vv, a.X_);
  int768 vvvR = int768_sub4(vvv, R);
  vvvR = int768_sub4(vvvR, R);
  int768 A = int768_sub4(int768_mul4(uu, a.Z_), vvvR);
  int768 X3 = int768_mul4(v, A);
  int768 vvvY1Z2 = int768_mul4(vvv, a.Y_);
  int768 Y3 = int768_sub4(int768_mul4(u, int768_sub4(R, A)), vvvY1Z2); 
  int768 Z3 = int768_mul4(vvv, a.Z_);
  
  res.X_ = X3;
  res.Y_ = Y3;
  res.Z_ = Z3;
  return res;
}

MNT_G1 G1_add6(MNT_G1 a, MNT_G1 b) {
  
}

MNT_G1 G1_mixed_add6(MNT_G1 a, MNT_G1 b) {
  
}



__kernel void mnt4753_fft(
    __global int768* x,
    __global int768* y,
    __global int768* pq,
    __global int768* omegas,
    __local int768* u,
    const unsigned int n,
    const unsigned int lgp,
    const unsigned int deg,
    const unsigned int max_deg) // 1=>radix2, 2=>radix4, 3=>radix8, ...
{
  uint32 lid = get_local_id(0);
  uint32 lsize = get_local_size(0);
  uint32 index = get_group_id(0);
  uint32 t = n >> deg;
  uint32 p = 1 << lgp;
  uint32 k = index & (p - 1);

  x += index;
  y += ((index - k) << deg) + k;

  uint32 count = 1 << deg; // 2^deg
  uint32 counth = count >> 1; // Half of count

  uint32 counts = count / lsize * lid;
  uint32 counte = counts + count / lsize;

  //////// ~30% of total time
  int768 twiddle = int768_pow_cached(omegas, (n >> lgp >> deg) * k);
  ////////

  //////// ~35% of total time
  int768 tmp = int768_pow(twiddle, counts);
  for(uint32 i = counts; i < counte; i++) {
    u[i] = int768_mul4(tmp, x[i*t]);
    tmp = int768_mul4(tmp, twiddle);
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  ////////

  //////// ~35% of total time
  uint32 pqshift = max_deg - deg;
  for(uint32 rnd = 0; rnd < deg; rnd++) {
    uint32 bit = counth >> rnd;
    for(uint32 i = counts >> 1; i < counte >> 1; i++) {
      uint32 di = i & (bit - 1);
      uint32 i0 = (i << 1) - di;
      uint32 i1 = i0 + bit;
      tmp = u[i0];
      u[i0] = int768_add4(u[i0], u[i1]);
      u[i1] = int768_sub4(tmp, u[i1]);
      if(di != 0) u[i1] = int768_mul4(pq[di << rnd << pqshift], u[i1]);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
  }
  ////////

  for(uint32 i = counts >> 1; i < counte >> 1; i++) {
    y[i*p] = u[bitreverse(i, deg)];
    y[(i+counth)*p] = u[bitreverse(i + counth, deg)];
  }
}

// Multi_exp functions
//

#define NUM_WORKS (224)
#define NUM_WINDOWS (96)
#define WINDOW_SIZE (8)
#define BUCKET_LEN ((1 << WINDOW_SIZE) - 1)

__kernel void G1_batched_lookup_multiexp(
    __global MNT_G1 *bases,
    __global MNT_G1 *buckets,
    __global MNT_G1 *results,
    __global int768 *exps,
    uint n) {

  uint32 gid = get_global_id(0);

  //bases += skip;
  buckets += BUCKET_LEN * gid;
  for(uint i = 0; i < BUCKET_LEN; i++) buckets[i] = G1_ZERO;

  uint len = (uint)ceil(n / (float)NUM_WORKS);
  uint32 nstart = len * (gid / NUM_WINDOWS);
  uint32 nend = min(nstart + len, n);

  uint bits = (gid % NUM_WINDOWS) * WINDOW_SIZE;
  //printf("%u\n",bits);
  ushort w = min((ushort)WINDOW_SIZE, (ushort)(768 - bits));

  MNT_G1 res = G1_ZERO;
  for(uint i = nstart; i < nend; i++) {
    uint ind = EXPONENT_get_bits(exps[i], bits, w);
    if(bits == 0 && ind == 1) res = G1_add4(res, bases[i]);
    else if(ind--) buckets[ind] = G1_add4(buckets[ind], bases[i]);
  }

  MNT_G1 acc = G1_ZERO;
  for(int j = BUCKET_LEN - 1; j >= 0; j--) {
    acc = G1_add4(acc, buckets[j]);
    res = G1_add4(res, acc);
  }

  results[gid] = res;
  //printf("%u\n", sizeof(MNT_G1));
}
