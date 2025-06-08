/* Implementation - Group F: Quaternion Operations */

/* Forward declarations for functions that may call each other */
FIXMATH_FUNC_ATTRS fix16_t qf16_dot(const qf16 *q, const qf16 *r);
FIXMATH_FUNC_ATTRS fix16_t qf16_norm(const qf16 *q);

/*---------------------------------------------------------------------------*/
/* BASIC QUATERNION OPERATIONS                                               */
/*---------------------------------------------------------------------------*/

/* Quaternion conjugate */
FIXMATH_FUNC_ATTRS void qf16_conj(qf16 *dest, const qf16 *q)
{
    dest->a = q->a;
    dest->b = -q->b;
    dest->c = -q->c;
    dest->d = -q->d;
}

/* Quaternion addition */
FIXMATH_FUNC_ATTRS void qf16_add(qf16 *dest, const qf16 *q, const qf16 *r)
{
    dest->a = q->a + r->a;
    dest->b = q->b + r->b;
    dest->c = q->c + r->c;
    dest->d = q->d + r->d;
}

/* Quaternion multiplication by scalar */
FIXMATH_FUNC_ATTRS void qf16_mul_s(qf16 *dest, const qf16 *q, fix16_t s)
{
    dest->a = fix16_mul(q->a, s);
    dest->b = fix16_mul(q->b, s);
    dest->c = fix16_mul(q->c, s);
    dest->d = fix16_mul(q->d, s);
}

/* Quaternion division by scalar */
FIXMATH_FUNC_ATTRS void qf16_div_s(qf16 *dest, const qf16 *q, fix16_t s)
{
    dest->a = fix16_div(q->a, s);
    dest->b = fix16_div(q->b, s);
    dest->c = fix16_div(q->c, s);
    dest->d = fix16_div(q->d, s);
}

/* Quaternion dot product */
FIXMATH_FUNC_ATTRS fix16_t qf16_dot(const qf16 *q, const qf16 *r)
{
    return fa16_dot(&q->a, &q->b - &q->a, &r->a, &r->b - &r->a, 4);    
}

/* Quaternion norm */
FIXMATH_FUNC_ATTRS fix16_t qf16_norm(const qf16 *q)
{
    return fa16_norm(&q->a, &q->b - &q->a, 4);
}

/* Quaternion normalization */
FIXMATH_FUNC_ATTRS void qf16_normalize(qf16 *dest, const qf16 *q)
{
    qf16_div_s(dest, q, qf16_norm(q));
}

/* Quaternion multiplication */
FIXMATH_FUNC_ATTRS void qf16_mul(qf16 *dest, const qf16 *q, const qf16 *r)
{
    qf16 tmp;
    fa16_unalias(dest, (void**)&q, (void**)&r, &tmp, sizeof(tmp));
    
    dest->a = fix16_mul(q->a, r->a) - fix16_mul(q->b, r->b) - fix16_mul(q->c, r->c) - fix16_mul(q->d, r->d);
    dest->b = fix16_mul(q->a, r->b) + fix16_mul(q->b, r->a) + fix16_mul(q->c, r->d) - fix16_mul(q->d, r->c);
    dest->c = fix16_mul(q->a, r->c) - fix16_mul(q->b, r->d) + fix16_mul(q->c, r->a) + fix16_mul(q->d, r->b);
    dest->d = fix16_mul(q->a, r->d) + fix16_mul(q->b, r->c) - fix16_mul(q->c, r->b) + fix16_mul(q->d, r->a);
}

/* Quaternion power */
FIXMATH_FUNC_ATTRS void qf16_pow(qf16 *dest, const qf16 *q, fix16_t power)
{
    fix16_t old_half_angle = fix16_acos(q->a);
    fix16_t new_half_angle = fix16_mul(old_half_angle, power);
    fix16_t multiplier = 0;
    
    if (old_half_angle > 10) // Guard against almost-zero divider
    {
        multiplier = fix16_div(fix16_sin(new_half_angle),
                               fix16_sin(old_half_angle));
    }
    
    dest->a = fix16_cos(new_half_angle);
    dest->b = fix16_mul(q->b, multiplier);
    dest->c = fix16_mul(q->c, multiplier);
    dest->d = fix16_mul(q->d, multiplier);
}

/* Quaternion weighted average */
FIXMATH_FUNC_ATTRS void qf16_avg(qf16 *dest, const qf16 *q1, const qf16 *q2, fix16_t weight)
{
    // z = sqrt((w1 - w2)^2 + 4 w1 w2 (q1' q2)^2
    // <=>
    // z = sqrt((2 w1 - 1)^2 + 4 w1 (1 - w1) (q1' q2)^2)
    fix16_t dot = qf16_dot(q1, q2);
    fix16_t z = fix16_sq(2 * weight - fix16_one)
            + fix16_mul(4 * weight, fix16_mul((fix16_one - weight), fix16_sq(dot)));
    z = fix16_sqrt(z);
    
    // q = 2 * w1 * (q1' q2) q1 + (w2 - w1 + z) q2
    // <=>
    // q = 2 * w1 * (q1' q2) q1 + (1 - 2 * w1 + z) q2
    qf16 tmp1;
    qf16_mul_s(&tmp1, q1, fix16_mul(2 * weight, dot));
    
    qf16 tmp2;
    qf16_mul_s(&tmp2, q2, fix16_one - 2 * weight + z);
    
    qf16_add(dest, &tmp1, &tmp2);
    qf16_normalize(dest, dest);
}

/*---------------------------------------------------------------------------*/
/* QUATERNION CONVERSION FUNCTIONS                                           */
/*---------------------------------------------------------------------------*/

/* Create unit quaternion from axis and angle */
FIXMATH_FUNC_ATTRS void qf16_from_axis_angle(qf16 *dest, const v3d *axis, fix16_t angle)
{
    angle /= 2;
    fix16_t scale = fix16_sin(angle);
    
    dest->a = fix16_cos(angle);
    dest->b = fix16_mul(axis->x, scale);
    dest->c = fix16_mul(axis->y, scale);
    dest->d = fix16_mul(axis->z, scale);
}

/* Convert unit quaternion to rotation matrix */
FIXMATH_FUNC_ATTRS void qf16_to_matrix(mf16 *dest, const qf16 *q)
{
    dest->rows = dest->columns = 3;
    dest->errors = 0;
    dest->data[0][0] = fix16_one - 2 * (fix16_sq(q->c) + fix16_sq(q->d));
    dest->data[1][1] = fix16_one - 2 * (fix16_sq(q->b) + fix16_sq(q->d));
    dest->data[2][2] = fix16_one - 2 * (fix16_sq(q->b) + fix16_sq(q->c));
    
    dest->data[1][0] = 2 * (fix16_mul(q->b, q->c) + fix16_mul(q->a, q->d));
    dest->data[0][1] = 2 * (fix16_mul(q->b, q->c) - fix16_mul(q->a, q->d));
    
    dest->data[2][0] = 2 * (fix16_mul(q->b, q->d) - fix16_mul(q->a, q->c));
    dest->data[0][2] = 2 * (fix16_mul(q->b, q->d) + fix16_mul(q->a, q->c));
    
    dest->data[2][1] = 2 * (fix16_mul(q->c, q->d) + fix16_mul(q->a, q->b));
    dest->data[1][2] = 2 * (fix16_mul(q->c, q->d) - fix16_mul(q->a, q->b));
} 