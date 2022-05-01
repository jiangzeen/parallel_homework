#include <iostream>
#include <chrono>
#include "openssl/ec.h"
#include <thread>
#include <vector>
#include "gmp.h"
#include "bigint.hpp"
using namespace std;
using namespace libff;
long long get_nsec_time();
void print_ECP(EC_GROUP *group, EC_POINT* p, BN_CTX *ctx);
void dbl_ECP(EC_GROUP *group, EC_POINT*P, BN_CTX *ctx);
void add_ECP(EC_GROUP *group, EC_POINT*P, BN_CTX *ctx);
void test_main(EC_GROUP *group, EC_POINT * Q, const vector<long>& res, BN_CTX *ctx);
//Sm2指定的参数，y2 = x3 + ax + b
#define _p  "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF"
#define _a  "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC"
#define _b  "28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93"
#define _n  "FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123"
#define _Gx "32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7"
#define _Gy "BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0"
// 计时器
typedef std::chrono::high_resolution_clock clock_type;//标准时钟：steady_clock。高分辨率时钟：high_resolution_clock
extern std::chrono::time_point<clock_type> TM_start{}, TM_end{};
#define Clock() clock_type::now()
#define Time(t_start,t_end) std::chrono::duration<double, ratio<1, 1>>(t_end - t_start).count() //时间差，类型 double, s
#define Timer(code) TM_start = Clock(); code; TM_end = Clock(); std::cout << Time(TM_start,TM_end) << " s\n"; //对code部分计时
// 循环测试
#define Loop(loop, code) Timer(for(long long i=0;i<loop;i++) {code;})

vector<long> naf_res;
vector<EC_POINT*> dbl_res;
EC_POINT * k_res;
int dbl_end;
template<mp_size_t n>
vector<long> find_naf(bigint<n> &scalar, const size_t window_size);
int main(int argc, char* argv[]) {
    BN_CTX  *ctx = BN_CTX_new();
    BIGNUM *p = BN_new(), *a = BN_new(), *b = BN_new(),
            *n = BN_new(), *Gx = BN_new(), *Gy = BN_new(), *h = BN_new();
    BN_hex2bn(&p, _p);
    BN_hex2bn(&a, _a);
    BN_hex2bn(&b, _b);
    BN_hex2bn(&n, _n);
    BN_hex2bn(&Gx, _Gx);
    BN_hex2bn(&Gy, _Gy);
    BN_set_word(h, 1);
    EC_GROUP *group = EC_GROUP_new(EC_GFp_mont_method());
    //设置曲线
    EC_GROUP_set_curve_GFp(group, p, a, b, ctx);
    //基点G
    EC_POINT * G = EC_POINT_new(group);
    EC_POINT_set_affine_coordinates_GFp(group, G, Gx, Gy, ctx);
    //设置基点
    EC_GROUP_set_generator(group, G, n, h);
    bigint<8> scalar1 = "8979845778459787512111878456789779879874512945877987456156464";
    EC_POINT *Q  = EC_POINT_new(group);
    EC_POINT_set_affine_coordinates(group, Q, BN_new(), BN_new(), ctx);
    print_ECP(group, Q, ctx);
    naf_res = find_naf(scalar1, 1);
    dbl_res = vector<EC_POINT*>(naf_res.size());

    long start = get_nsec_time();
    //
    thread process1(dbl_ECP, group, G, ctx);
    thread process2(add_ECP, group, G, ctx);
    process1.join();
    process2.join();
    return 0;
}

//
void dbl_ECP(EC_GROUP *group, EC_POINT*P, BN_CTX *ctx){
    EC_POINT * G = EC_POINT_new(group);
    EC_POINT_copy(G, P);
    for (int i = 0; i < naf_res.size(); ++i) {
        if (naf_res[i] == 1) {
            EC_POINT * newP = EC_POINT_new(group);
            EC_POINT_copy(newP, G);
            dbl_res[i] = newP;
        } else if (naf_res[i] == -1) {
            EC_POINT * newP = EC_POINT_new(group);
            EC_POINT_copy(newP, G);
            EC_POINT_invert(group, newP, ctx);
            dbl_res[i] = newP;
        }
        EC_POINT_dbl(group, G, G, ctx);
    }
    dbl_end = 1;
}
void add_ECP(EC_GROUP *group, EC_POINT*P, BN_CTX *ctx){
    EC_POINT * G = EC_POINT_new(group);
    EC_POINT_set_affine_coordinates(group, G, BN_new(), BN_new(), ctx);
    // EC_POINT_copy(G, P);
    int i = 0;
    while (dbl_end != 1) {
        while (dbl_res[i] == nullptr) {
            // loop while no res
            if (dbl_res[i] != nullptr) break;
        }
        EC_POINT_add(group, G, dbl_res[i], G, ctx);
        i++;
    }
    while (dbl_res[i] != nullptr) {
        EC_POINT_add(group, G, dbl_res[i], G, ctx);
        i++;
    }
    k_res = EC_POINT_new(group);
    EC_POINT_copy(k_res, G);
}



// print point
void print_ECP(EC_GROUP *group, EC_POINT* p, BN_CTX *ctx)
{
    printf("struct EC_POINT\n");

    BIGNUM *x = BN_new(), *y = BN_new();

    //s = 04 || x || y
    /*char* s = EC_POINT_point2hex(group, p, POINT_CONVERSION_UNCOMPRESSED, ctx);
    printf("\tPoint = %s\n", s);*/

    EC_POINT_get_affine_coordinates_GFp(group, p, x, y, ctx);

    char *s1 = BN_bn2hex(x);
    char *s2 = BN_bn2hex(y);
    printf("\tx = %s\n", s1);
    printf("\ty = %s\n", s2);
    BN_free(x);
    BN_free(y);
}

// calc naf form
template<mp_size_t n>
vector<long> find_naf(bigint<n> &scalar, const size_t window_size) {
    const size_t length = scalar.max_bits(); // upper bound
    std::vector<long> res(length+1);
    bigint<n> c = scalar;
    long j = 0;
    while (!c.is_zero())
    {
        long u;
        if ((c.data[0] & 1) == 1)
        {
            u = c.data[0] % (1u << (window_size+1));
            if (u > (1 << window_size))
            {
                u = u - (1 << (window_size+1));
            }

            if (u > 0)
            {
                mpn_sub_1(c.data, c.data, n, u);
            }
            else
            {
                mpn_add_1(c.data, c.data, n, -u);
            }
        }
        else
        {
            u = 0;
        }
        res[j] = u;
        ++j;

        mpn_rshift(c.data, c.data, n, 1); // c = c/2
    }
    return res;
}

long long get_nsec_time()
{
    auto timepoint = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(timepoint.time_since_epoch()).count();
}

void test_main(EC_GROUP *group, EC_POINT * Q, const vector<long>& res, BN_CTX *ctx) {
    EC_POINT * M = EC_POINT_new(group);
    EC_POINT_copy(M, Q);
    for (long i = res.size() - 1; i >= 0; i--) {
        EC_POINT_dbl(group, Q, Q, ctx);
        if (res[i] == 1) {
            EC_POINT_add(group, Q, Q, M, ctx);
        } else if (res[i] == -1) {
            EC_POINT_invert(group, M, ctx);
            EC_POINT_add(group, Q, Q, M, ctx);
            EC_POINT_invert(group, M, ctx);
        }
    }
}