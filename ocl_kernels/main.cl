// FinalityLabs - 2019
// fixed 768 size prime-field arithmetic library (add, sub, mul, pow)
// Montgomery reduction parameters:
// B = 2^32 (Because our digits are uint32)


typedef uint uint32;
typedef ulong uint64;

typedef uint32 limb;
typedef uint64 limb2;

typedef struct {
  limb v[24];
} int768;

#define FIELD_LIMBS (24)
#define LIMB_BITS (32)
#define LIMB_MAX (0xffffffff)

// Montgomery form of 1 = (1 * R mod P)
//
#define mnt4753_ONE ((int768){{0xd9dc6f42,0x98a8ecab,0x5a034686,0x91cd31c6,0xcd14572e,0x97c3e4a0,0xc788b601,0x79589819,0x2108976f,0xed269c94,0xcf031d68,0x1e0f4d8a,0x13338559,0x320c3bb7,0xd2f00a62,0x598b4302,0xfd8ca621,0x4074c9cb,0x3865e88c,0xfa47edb,0x1ff9a195,0x95455fb3,0x9ec8e242,0x7b47}})
//
//#define mnt6753_ONE ((int768){0x,0x,0x,0x,0x,0x,0x,0x,0x,0x,0x,0x})

#define mnt4753_ZERO ((int768){{0}})
#define mnt6753_ZERO ((int768){{0}})

//#define mnt4753_INV_Fr ((ulong)0xc90776e23fffffff)
#define mnt4753_INV_Fq ((uint){{0xe45e7fff}})
//#define mnt6753_INV_Fr ((ulong)0xf2044cfbe45e7fff)
#define mnt6753_INV_Fq ((uint){{0x3fffffff}})

#define mnt4753_Q ((int768){{0x245e8001,0x5e9063de,0x2cdd119f,0xe39d5452,0x9ac425f0,0x63881071,0x767254a4,0x685acce9,0xcb537e38,0xb80f0da5,0xf218059d,0xb117e776,0xa15af79d,0x99d124d9,0xe8a0ed8d,0x07fdb925,0x6c97d873,0x5eb7e8f9,0x5b8fafed,0xb7f99750,0xeee2cdad,0x10229022,0x2d92c411,0x1c4c6}})
#define mnt6753_Q ((int768){{0x40000001,0xd90776e2,0xfa13a4f,0x4ea09917,0x3f005797,0xd6c381bc,0x34993aa4,0xb9dff976,0x29212636,0x3eebca94,0xc859a99b,0xb26c5c28,0xa15af79d,0x99d124d9,0xe8a0ed8d,0x7fdb925,0x6c97d873,0x5eb7e8f9,0x5b8fafed,0xb7f99750,0xeee2cdad,0x10229022,0x2d92c411,0x1c4c6}})

limb bitreverse(limb n, limb bits) {
  limb r = 0;
  for(int i = 0; i < bits; i++) {
    r = (r << 1) | (n & 1);
    n >>= 1;
  }
  return r;
}

void print(int768 v) {
  printf("%u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u\n",
    v.v[23],v.v[22],v.v[21],v.v[20],v.v[19],v.v[18],v.v[17],v.v[16],v.v[15],v.v[14],v.v[13],v.v[12], v.v[11],v.v[10],v.v[9],v.v[8],v.v[7],v.v[6],v.v[5],v.v[4],v.v[3],v.v[2],v.v[1],v.v[0]);
}


bool get_bit(int768 l, uint i) {

  if(i < 32)
    return (l.v[23] >> (31 - i)) & 1;
  else if(i < 64)
    return (l.v[22] >> (63 - i)) & 1;
  else if(i < 96)
    return (l.v[21] >> (95 - i)) & 1;
  else if(i < 128)
    return (l.v[20] >> (127 - i)) & 1;
  else if(i < 160)
    return (l.v[19] >> (159 - i)) & 1;
  else if(i < 192)
    return (l.v[18] >> (192 - i)) & 1;
  else if(i < 224)
    return (l.v[17] >> (223 - i)) & 1;
  else if(i < 256)
    return (l.v[16] >> (255 - i)) & 1;
  else if(i < 288)
    return (l.v[15] >> (287 - i)) & 1;
  else if(i < 320)
    return (l.v[14] >> (319 - i)) & 1;
  else if(i < 352)
    return (l.v[13] >> (351 - i)) & 1;
  else if(i < 384)
    return (l.v[12] >> (383 - i)) & 1;
  else if(i < 416)
    return (l.v[11] >> (415 - i)) & 1;
  else if(i < 448)
    return (l.v[10] >> (447 - i)) & 1;
  else if(i < 480)
    return (l.v[9] >> (479 - i)) & 1;
  else if(i < 512)
    return (l.v[8] >> (511 - i)) & 1;
  else if(i < 544)
    return (l.v[7] >> (543 - i)) & 1;
  else if(i < 576)
    return (l.v[6] >> (575 - i)) & 1;
  else if(i < 608)
    return (l.v[5] >> (607 - i)) & 1;
  else if(i < 640)
    return (l.v[4] >> (639 - i)) & 1;
  else if(i < 672)
    return (l.v[3] >> (671 - i)) & 1;
  else if(i < 704)
    return (l.v[2] >> (703 - i)) & 1;
  else if(i < 736)
    return (l.v[1] >> (735 - i)) & 1;
  else
    return (l.v[0] >> (767 - i)) & 1;
}

uint get_bits(int768 l, uint skip, uint window) {
  uint ret = 0;
  for(uint i = 0; i < window; i++) {
    ret <<= 1;
    ret |= get_bit(l, skip + i);
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

// Adds `num` to `i`th digit of `res` and propagates carry in case of overflow
void add_digit(limb *res, limb num) {
  limb old = *res;
  *res += num;
  if(*res < old) {
    res++;
    while(++(*(res++)) == 0);
  }
}

limb mac_with_carry(limb a, limb b, limb c, limb *carry) {
  limb lo = a * b;
  limb hi = mul_hi(a, b);
  hi += lo + c < lo; lo += c;
  hi += lo + *carry < lo; lo += *carry;
  *carry = hi;
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
    limb2 sum = (limb2)a.v[i] + b.v[i] + carry;
    a.v[i] = sum & LIMB_MAX;
    carry = sum >> LIMB_BITS;
    //"this implementation fails for some inputs"
    //ulong old = a.v[i];
    //a.v[i] += b.v[i] + carry;
    //carry = carry ? old >= a.v[i] : old > a.v[i];
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
    ulong old = a.v[i];
    a.v[i] -= b.v[i] + borrow;
    borrow = borrow ? old <= a.v[i] : old < a.v[i];
  }
  return a;
}

int768 int768_reduce(ulong *limbs) {
  // Montgomery reduction
  bool carry2 = 0;
  for(uchar i = 0; i < FIELD_LIMBS; i++) {
    limb u = mnt4753_INV_Fq * limbs[i];
    limb carry = 0;
    for(uchar j = 0; j < FIELD_LIMBS; j++) {
      limb2 product = (limb2)u * mnt4753_Q.v[j] + limbs[i + j] + carry;
      limbs[i + j] = product & LIMB_MAX;
      carry = product >> LIMB_BITS;
    }
    limb2 sum = (limb2)limbs[i + FIELD_LIMBS] + carry + carry2;
    limbs[i + FIELD_LIMBS] = sum & LIMB_MAX;
    carry2 = sum >> LIMB_BITS;
  }

  // Divide by R
  int768 result;
  // this breaks amd compiler
  for(uchar i = 0; i < FIELD_LIMBS; i++) result.v[i] = limbs[i+FIELD_LIMBS];

  if(int768_gte(result, mnt6753_Q))
    result = int768_sub_(result, mnt6753_Q);

  return result;
}

// Modular multiplication
int768 int768_mul4(int768 a, int768 b) {
  // long multiplication
  limb res[FIELD_LIMBS * 2] = {0};
  for(uint32 i = 0; i < FIELD_LIMBS; i++) {
    limb carry = 0;
    for(uint32 j = 0; j < FIELD_LIMBS; j++) {
      limb2 product = (limb2)a.v[i] * b.v[j] + res[i + j] + carry;
      res[i + j] = product & LIMB_MAX;
      carry = product >> LIMB_BITS;
      //res[i + j] = mac_with_carry(a.v[i], b.v[j], res[i + j], &carry);
    }
    res[i + FIELD_LIMBS] = carry;
  }
  // this breaks amd compiler
  //return int768_reduce(res);

  bool carry2 = 0;
  for(uchar i = 0; i < FIELD_LIMBS; i++) {
    limb u = mnt4753_INV_Fq * res[i];
    limb carry = 0;
    for(uchar j = 0; j < FIELD_LIMBS; j++) {
      limb2 product = (limb2)u * mnt4753_Q.v[j] + res[i + j] + carry;
      res[i + j] = product & LIMB_MAX;
      carry = product >> LIMB_BITS;
    }
    limb2 sum = (limb2)res[i + FIELD_LIMBS] + carry + carry2;
    res[i + FIELD_LIMBS] = sum & LIMB_MAX;
    carry2 = sum >> LIMB_BITS;
  }

  // Divide by R
  int768 result;
  for(uchar i = 0; i < FIELD_LIMBS; i++) result.v[i] = res[i+FIELD_LIMBS];

  if(int768_gte(result, mnt4753_Q))
    result = int768_sub_(result, mnt4753_Q);

  return result;
}

// Modular multiplication
int768 int768_mul6(int768 a, int768 b) {
  // long multiplication
  limb res[FIELD_LIMBS * 2] = {0};
  for(uint32 i = 0; i < FIELD_LIMBS; i++) {
    limb carry = 0;
    for(uint32 j = 0; j < FIELD_LIMBS; j++) {
      limb2 product = (limb2)a.v[i] * b.v[j] + res[i + j] + carry;
      res[i + j] = product & LIMB_MAX;
      carry = product >> LIMB_BITS;
      //res[i + j] = mac_with_carry(a.v[i], b.v[j], res[i + j], &carry);
    }
    res[i + FIELD_LIMBS] = carry;
  }
  // this breaks amd compiler
  //return int768_reduce(res);

  bool carry2 = 0;
  for(uchar i = 0; i < FIELD_LIMBS; i++) {
    limb u = mnt6753_INV_Fq * res[i];
    limb carry = 0;
    for(uchar j = 0; j < FIELD_LIMBS; j++) {
      limb2 product = (limb2)u * mnt6753_Q.v[j] + res[i + j] + carry;
      res[i + j] = product & LIMB_MAX;
      carry = product >> LIMB_BITS;
    }
    limb2 sum = (limb2)res[i + FIELD_LIMBS] + carry + carry2;
    res[i + FIELD_LIMBS] = sum & LIMB_MAX;
    carry2 = sum >> LIMB_BITS;
  }

  // Divide by R
  int768 result;
  for(uchar i = 0; i < FIELD_LIMBS; i++) result.v[i] = res[i+FIELD_LIMBS];

  if(int768_gte(result, mnt6753_Q))
    result = int768_sub_(result, mnt6753_Q);

  return result;
}

// Modular negation
int768 int768_neg4(int768 a) {
  return int768_sub_(mnt4753_Q, a);
}

int768 int768_neg6(int768 a) {
  return int768_sub_(mnt6753_Q, a);
}

// Modular subtraction
int768 int768_sub4(int768 a, int768 b) {
  int768 res = int768_sub_(a, b);
  if(!int768_gte(a, b)) res = int768_add_(res, mnt4753_Q);
  return res;
}

int768 int768_sub6(int768 a, int768 b) {
  int768 res = int768_sub_(a, b);
  if(!int768_gte(a, b)) res = int768_add_(res, mnt6753_Q);
  return res;
}

// Modular addition
int768 int768_add4(int768 a, int768 b) {
  //return int768_sub4(a, int768_neg4(b));
  int768 tmp = int768_neg4(b);
  int768 res = int768_sub_(a, tmp);
  if(!int768_gte(a, tmp)) res = int768_add_(res, mnt4753_Q);
  return res;
}

int768 int768_add6(int768 a, int768 b) {
  //return int768_sub4(a, int768_neg6(b));
  int768 tmp = int768_neg6(b);
  int768 res = int768_sub_(a, tmp);
  if(!int768_gte(a, tmp)) res = int768_add_(res, mnt6753_Q);
  return res;
}

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
#define non_residue ((int768){{0xa3162657,0xa4e2d91f,0xb935ff7,0xbc938a1c,0x99bbfb8a,0x8a5a6ad5,0xbe9a4027,0xf06f5292,0x4b7535ff,0xe2c8ca94,0xace06d7a,0x737f39a7,0x158cdead,0xbd2b99bf,0xfc4dbe53,0x74193bb2,0x9a5ce658,0x29c6846f,0xca7dbf57,0xa36dab30,0xd304cb88,0x641e2baf,0x877b312e,0xf450}});
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
  //print(residue);
  int768 A = _b.c0;
  int768 a = _a.c0;
  int768 B = _b.c1;
  int768 b = _a.c1;

  //  const my_Fp
  //      &A = other.c0, &B = other.c1,
  //      &a = this->c0, &b = this->c1;
  //  const my_Fp aA = a * A;
  //  const my_Fp bB = b * B;


  int768 aA = int768_mul4(a, A);
  int768 bB = int768_mul4(b, B);

  Fq2 res = Fq2_ZERO;

  //res.c0 = int768_add4(aa, int768_mul4(residue, bb));
  //res.c0 = int768_add4(aA, int768_mul4(residue, bB));
  res.c0 = int768_add4(int768_mul4(_a.c0, _b.c0), int768_mul4(residue, int768_mul4(_a.c1, _b.c1)));
  
  
  // Sub(Sub(Mul(Add(x.a, x.b), Add(y.a, y.b)), A), B)
  //return Fp2_model<n,modulus>(aA + non_residue * bB,
  //                            (a + b)*(A+B) - aA - bB);

  //int768 v4 = int768_add4(a, b);
  //int768 v3 = int768_add4(A, B);
  //int768 v2 = int768_mul4(v4, v3);
  //int768 v1 = int768_sub4(v2, aA);
  //int768 v0 = int768_sub4(v1, bB);
  
  //res.c1 = v0;
  res.c1 = int768_mul4(int768_add4(_a.c0, _a.c1), int768_add4(_b.c0, _b.c1));
  res.c1 = int768_sub4(res.c1, aA);
  res.c1 = int768_sub4(res.c1, bB);
  //res.c1 = int768_sub4(res.c1, int768_mul4(a, A));
  return res;
}


// Fq3 arithmetics
//

typedef struct {
  int768 c0;
  int768 c1;
  int768 c2;
} Fq3;

#define Fq3_ZERO ((Fq3){mnt6753_ZERO, mnt6753_ZERO, mnt6753_ZERO})
#define Fq3_ONE ((Fq3){mnt6753_ONE, mnt6753_ZERO, mnt6753_ZERO})

// Montgomery non_residue
#define non_residue ((int768){{0xfff9c7d4,0x4768931c,0xada96ca0,0xc45e46d6,0xb3c0107,0x479b0bdb,0x10f8d41b,0x362a0896,0xc8a91aaf,0xdbafcec2,0xf9d96a06,0x78428b0f,0x9080c353,0xf2e4472a,0x3f0e971c,0xc9006ed3,0xbdb7288,0x794d9d1,0xb5419e2c,0x3c1e44ca,0x81f4560c,0x49b5fc6c,0x777c30ba,0x1c287}});
// Integer non_residue
//#define non_residue ((int768){{0xD}});


Fq3 int768_fq3_add(Fq3 x, Fq3 y) {
  Fq3 res = Fq3_ZERO;
  res.c0 = int768_add6(x.c0, y.c0); 
  res.c1 = int768_add6(x.c1, y.c1);
  res.c2 = int768_add6(x.c2, y.c2);
  return res;
}

Fq3 int768_fq3_mul(Fq3 a, Fq3 b) {
  Fq3 res = Fq3_ZERO;
  int768 residue = non_residue;

  int768 a0_b0 = int768_mul6(a.c0, b.c0);
  int768 a0_b1 = int768_mul6(a.c0, b.c1); 
  int768 a0_b2 = int768_mul6(a.c0, b.c2);

  int768 a1_b0 = int768_mul6(a.c1, b.c0);
  int768 a1_b1 = int768_mul6(a.c1, b.c1); 
  int768 a1_b2 = int768_mul6(a.c1, b.c2);

  int768 a2_b0 = int768_mul6(a.c2, b.c0);
  int768 a2_b1 = int768_mul6(a.c2, b.c1); 
  int768 a2_b2 = int768_mul6(a.c2, b.c2);

  res.c0 = int768_add6(a0_b0, int768_mul6(residue, int768_add6(a1_b2, a2_b1)));
  res.c1 = int768_add6(a0_b1, int768_add6(a1_b0, int768_mul6(residue, a2_b2)));
  res.c2 = int768_add6(a0_b2, int768_add6(a1_b1, a2_b0));
  return res;
}

// eliptic curve arithmetics
//

typedef struct {
  int768 X_;
  int768 Y_;
  int768 Z_;
} MNT_G1;

// affine coord zero
#define G1_ZERO ((MNT_G1){mnt4753_ZERO, mnt4753_ONE, mnt4753_ZERO})
#define G1_COEFF_A ((int768){{0xb3b8de84,0x3151d957,0xb4068d0d,0x239a638c,0x9a28ae5d,0x2f87c941,0x8f116c03,0xf2b13033,0x42112ede,0xda4d3928,0x9e063ad1,0x3c1e9b15,0x26670ab2,0x6418776e,0xa5e014c4,0xb3168605,0xfb194c42,0x80e99397,0x70cbd118,0x1f48fdb6,0x3ff3432a,0x2a8abf66,0x3d91c485,0xf68f}})

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
    int768 w = int768_add4(int768_mul4(G1_COEFF_A, ZZ), TXX);
    int768 Y1_Z1 = int768_mul4(a.Y_, a.Z_);
    int768 s = int768_add4(Y1_Z1, Y1_Z1);
    int768 ss = int768_mul4(s, s);
    int768 sss = int768_mul4(s, ss);
    int768 R = int768_mul4(a.Y_, s);
    int768 RR = int768_mul4(R, R);
    int768 XRR = int768_add4(a.X_, RR);
    int768 B = int768_sub4(XRR, int768_sub4(XX, RR));
    int768 h = int768_sub4(int768_mul4(w, w), int768_add4(B, B));
    int768 X3 = int768_mul4(h, s);
    int768 Y3 = int768_mul4(w, int768_sub4(int768_sub4(B, h), int768_add4(R,R)));
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

  MNT_G1 res = G1_ZERO;

  int768 XX = int768_mul4(a.X_, a.X_); // todo special case squaring
  int768 ZZ = int768_mul4(a.Z_, a.Z_);
  int768 TXX = int768_add4(XX, XX);
  TXX = int768_add4(TXX, XX);
  int768 w = int768_add4(int768_mul4(G1_COEFF_A, ZZ), TXX);
  int768 Y1_Z1 = int768_mul4(a.Y_, a.Z_);
  int768 s = int768_add4(Y1_Z1, Y1_Z1);
  int768 ss = int768_mul4(s, s);
  int768 sss = int768_mul4(s, ss);
  int768 R = int768_mul4(a.Y_, s);
  int768 RR = int768_mul4(R, R);
  int768 XRR = int768_add4(a.X_, RR);
  int768 B = int768_sub4(XRR, int768_sub4(XX, RR));
  int768 h = int768_sub4(int768_mul4(w, w), int768_add4(B, B));
  int768 X3 = int768_mul4(h, s);
  int768 Y3 = int768_mul4(w, int768_sub4(int768_sub4(B, h), int768_add4(R,R)));
  res.X_ = X3;
  res.Y_ = Y3;
  res.Z_ = sss;
  return res;
}


MNT_G1 G1_mixed_add4(MNT_G1 a, MNT_G1 b) {
  
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

#define WINDOW_SIZE (4)
#define TABLE_SIZE ((1 << WINDOW_SIZE) - 1)

__kernel void G1_generate_table(
    __global MNT_G1 *bases,
    __global bool *dm,
    uint skip,
    uint n) {
  uint32 gid = get_global_id(0);
  bases += skip + gid;
  for(uint j = 1; j < TABLE_SIZE; j++)
    if(j & 1) bases[j * n] = G1_double4(bases[(j >> 1) * n]);
    else bases[j * n] = G1_add4(bases[(j - 1) * n], bases[0]);
}

__kernel void G1_batched_lookup_multiexp(
    __global MNT_G1 *bases,
    __global MNT_G1 *results,
    __global ulong4 *exps,
    __global bool *dm,
    uint skip,
    uint n) {
  uint32 work = get_global_id(0);
  uint32 works = get_global_size(0);

  uint len = (uint)ceil(n / (float)works);
  uint32 nstart = len * work;
  uint32 nend = min(nstart + len, n);

  bases += skip;

  MNT_G1 p = G1_ZERO;
  ushort bits = 0;
  // not sure if extra 15 bits matters
  while(bits < 768) {
    ushort w = min((ushort)WINDOW_SIZE, (ushort)(768 - bits));
    for(uint j = 0; j < w; j++)
      p = G1_double4(p);
    for(uint j = nstart; j < nend; j++) {
      uint ind = get_bits(exps[j], bits, w);
      //uint ind = 0;
      if(ind)
        p = G1_add4(p, bases[j + (ind - 1) * n]);
    }
    bits += w;
  }
  results[work] = p;
}
