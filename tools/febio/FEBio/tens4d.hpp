// NOTE: This file is automatically included from tens4d.h
// Users should not include this file manually!


// operator +
inline tens4ds tens4ds::operator + (const tens4ds& t) const
{
	tens4ds s;
//	for (int i=0; i<NNZ; i++)
//		s.d[i] = d[i] + t.d[i];
	s.d[0] = d[0] + t.d[0];
	s.d[1] = d[1] + t.d[1];
	s.d[2] = d[2] + t.d[2];
	s.d[3] = d[3] + t.d[3];
	s.d[4] = d[4] + t.d[4];
	s.d[5] = d[5] + t.d[5];
	s.d[6] = d[6] + t.d[6];
	s.d[7] = d[7] + t.d[7];
	s.d[8] = d[8] + t.d[8];
	s.d[9] = d[9] + t.d[9];
	s.d[10] = d[10] + t.d[10];
	s.d[11] = d[11] + t.d[11];
	s.d[12] = d[12] + t.d[12];
	s.d[13] = d[13] + t.d[13];
	s.d[14] = d[14] + t.d[14];
	s.d[15] = d[15] + t.d[15];
	s.d[16] = d[16] + t.d[16];
	s.d[17] = d[17] + t.d[17];
	s.d[18] = d[18] + t.d[18];
	s.d[19] = d[19] + t.d[19];
	s.d[20] = d[20] + t.d[20];
	return s;
}

// operator -
inline tens4ds tens4ds::operator - (const tens4ds& t) const
{
	tens4ds s;
//	for (int i=0; i<NNZ; i++)
//		s.d[i] = d[i] - t.d[i];
	s.d[0] = d[0] - t.d[0];
	s.d[1] = d[1] - t.d[1];
	s.d[2] = d[2] - t.d[2];
	s.d[3] = d[3] - t.d[3];
	s.d[4] = d[4] - t.d[4];
	s.d[5] = d[5] - t.d[5];
	s.d[6] = d[6] - t.d[6];
	s.d[7] = d[7] - t.d[7];
	s.d[8] = d[8] - t.d[8];
	s.d[9] = d[9] - t.d[9];
	s.d[10] = d[10] - t.d[10];
	s.d[11] = d[11] - t.d[11];
	s.d[12] = d[12] - t.d[12];
	s.d[13] = d[13] - t.d[13];
	s.d[14] = d[14] - t.d[14];
	s.d[15] = d[15] - t.d[15];
	s.d[16] = d[16] - t.d[16];
	s.d[17] = d[17] - t.d[17];
	s.d[18] = d[18] - t.d[18];
	s.d[19] = d[19] - t.d[19];
	s.d[20] = d[20] - t.d[20];
	return s;
}

// operator *
inline tens4ds tens4ds::operator * (double g) const
{
	tens4ds s;
//	for (int i=0; i<NNZ; i++)
//		s.d[i] = g*d[i];
	s.d[0] = g*d[0];
	s.d[1] = g*d[1];
	s.d[2] = g*d[2];
	s.d[3] = g*d[3];
	s.d[4] = g*d[4];
	s.d[5] = g*d[5];
	s.d[6] = g*d[6];
	s.d[7] = g*d[7];
	s.d[8] = g*d[8];
	s.d[9] = g*d[9];
	s.d[10] = g*d[10];
	s.d[11] = g*d[11];
	s.d[12] = g*d[12];
	s.d[13] = g*d[13];
	s.d[14] = g*d[14];
	s.d[15] = g*d[15];
	s.d[16] = g*d[16];
	s.d[17] = g*d[17];
	s.d[18] = g*d[18];
	s.d[19] = g*d[19];
	s.d[20] = g*d[20];
	return s;
}

// operator /
inline tens4ds tens4ds::operator / (double g) const
{
	tens4ds s;
//	for (int i=0; i<NNZ; i++)
//		s.d[i] = d[i]/g;
	s.d[0] = d[0]/g;
	s.d[1] = d[1]/g;
	s.d[2] = d[2]/g;
	s.d[3] = d[3]/g;
	s.d[4] = d[4]/g;
	s.d[5] = d[5]/g;
	s.d[6] = d[6]/g;
	s.d[7] = d[7]/g;
	s.d[8] = d[8]/g;
	s.d[9] = d[9]/g;
	s.d[10] = d[10]/g;
	s.d[11] = d[11]/g;
	s.d[12] = d[12]/g;
	s.d[13] = d[13]/g;
	s.d[14] = d[14]/g;
	s.d[15] = d[15]/g;
	s.d[16] = d[16]/g;
	s.d[17] = d[17]/g;
	s.d[18] = d[18]/g;
	s.d[19] = d[19]/g;
	s.d[20] = d[20]/g;
	return s;
}

// assignment operator +=
inline tens4ds& tens4ds::operator += (const tens4ds& t)
{
//	for (int i=0; i<NNZ; i++)
//		d[i] += t.d[i];
	d[0] += t.d[0];
	d[1] += t.d[1];
	d[2] += t.d[2];
	d[3] += t.d[3];
	d[4] += t.d[4];
	d[5] += t.d[5];
	d[6] += t.d[6];
	d[7] += t.d[7];
	d[8] += t.d[8];
	d[9] += t.d[9];
	d[10] += t.d[10];
	d[11] += t.d[11];
	d[12] += t.d[12];
	d[13] += t.d[13];
	d[14] += t.d[14];
	d[15] += t.d[15];
	d[16] += t.d[16];
	d[17] += t.d[17];
	d[18] += t.d[18];
	d[19] += t.d[19];
	d[20] += t.d[20];
	return (*this);
}

// assignment operator -=
inline tens4ds& tens4ds::operator -= (const tens4ds& t)
{
//	for (int i=0; i<NNZ; i++)
//		d[i] -= t.d[i];
	d[0] -= t.d[0];
	d[1] -= t.d[1];
	d[2] -= t.d[2];
	d[3] -= t.d[3];
	d[4] -= t.d[4];
	d[5] -= t.d[5];
	d[6] -= t.d[6];
	d[7] -= t.d[7];
	d[8] -= t.d[8];
	d[9] -= t.d[9];
	d[10] -= t.d[10];
	d[11] -= t.d[11];
	d[12] -= t.d[12];
	d[13] -= t.d[13];
	d[14] -= t.d[14];
	d[15] -= t.d[15];
	d[16] -= t.d[16];
	d[17] -= t.d[17];
	d[18] -= t.d[18];
	d[19] -= t.d[19];
	d[20] -= t.d[20];
	return (*this);
}

// assignment operator *=
inline tens4ds& tens4ds::operator *= (double g)
{
//	for (int i=0; i<NNZ; i++)
//		d[i] *= g;
	d[0] *= g;
	d[1] *= g;
	d[2] *= g;
	d[3] *= g;
	d[4] *= g;
	d[5] *= g;
	d[6] *= g;
	d[7] *= g;
	d[8] *= g;
	d[9] *= g;
	d[10] *= g;
	d[11] *= g;
	d[12] *= g;
	d[13] *= g;
	d[14] *= g;
	d[15] *= g;
	d[16] *= g;
	d[17] *= g;
	d[18] *= g;
	d[19] *= g;
	d[20] *= g;
	return (*this);
}

// assignment operator /=
inline tens4ds& tens4ds::operator /= (double g)
{
//	for (int i=0; i<NNZ; i++)
//		d[i] /= g;
	d[0] /= g;
	d[1] /= g;
	d[2] /= g;
	d[3] /= g;
	d[4] /= g;
	d[5] /= g;
	d[6] /= g;
	d[7] /= g;
	d[8] /= g;
	d[9] /= g;
	d[10] /= g;
	d[11] /= g;
	d[12] /= g;
	d[13] /= g;
	d[14] /= g;
	d[15] /= g;
	d[16] /= g;
	d[17] /= g;
	d[18] /= g;
	d[19] /= g;
	d[20] /= g;
	return (*this);
}

// intialize to zero
inline void tens4ds::zero()
{
	d[0] = d[1] = d[2] = d[3] = d[4] = d[5] = d[6] = d[7] = d[8] = d[9] =
	d[10] = d[11] = d[12] = d[13] = d[14] = d[15] = d[16] = d[17] = d[18] = d[19] = d[20] = 0;
}

// extract 6x6 matrix
inline void tens4ds::extract(double D[6][6])
{
	D[0][0] = d[0];  D[0][1] = d[1];  D[0][2] = d[3];  D[0][3] = d[6];  D[0][4] = d[10]; D[0][5] = d[15];
	D[1][0] = d[1];  D[1][1] = d[2];  D[1][2] = d[4];  D[1][3] = d[7];  D[1][4] = d[11]; D[1][5] = d[16];
	D[2][0] = d[3];  D[2][1] = d[4];  D[2][2] = d[5];  D[2][3] = d[8];  D[2][4] = d[12]; D[2][5] = d[17];
	D[3][0] = d[6];  D[3][1] = d[7];  D[3][2] = d[8];  D[3][3] = d[9];  D[3][4] = d[13]; D[3][5] = d[18];
	D[4][0] = d[10]; D[4][1] = d[11]; D[4][2] = d[12]; D[4][3] = d[13]; D[4][4] = d[14]; D[4][5] = d[19];
	D[5][0] = d[15]; D[5][1] = d[16]; D[5][2] = d[17]; D[5][3] = d[18]; D[5][4] = d[19]; D[5][5] = d[20];
}

//-----------------------------------------------------------------------------
// (a dyad1s a)_ijkl = a_ij a_kl
inline tens4ds dyad1s(const mat3ds& a)
{
	tens4ds c;
	c.d[ 0] = a.xx()*a.xx();
	c.d[ 1] = a.xx()*a.yy();
	c.d[ 3] = a.xx()*a.zz();
	c.d[ 6] = a.xx()*a.xy();
	c.d[10] = a.xx()*a.yz();
	c.d[15] = a.xx()*a.xz();

	c.d[ 2] = a.yy()*a.yy();
	c.d[ 4] = a.yy()*a.zz();
	c.d[ 7] = a.yy()*a.xy();
	c.d[11] = a.yy()*a.yz();
	c.d[16] = a.yy()*a.xz();

	c.d[ 5] = a.zz()*a.zz();
	c.d[ 8] = a.zz()*a.xy();
	c.d[12] = a.zz()*a.yz();
	c.d[17] = a.zz()*a.xz();

	c.d[ 9] = a.xy()*a.xy();
	c.d[13] = a.xy()*a.yz();
	c.d[18] = a.xy()*a.xz();

	c.d[14] = a.yz()*a.yz();
	c.d[19] = a.yz()*a.xz();

	c.d[20] = a.xz()*a.xz();
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad1s b)_ijkl = a_ij b_kl + b_ij a_kl
inline tens4ds dyad1s(const mat3ds& a, const mat3ds& b)
{
	tens4ds c;
	c.d[0] = 2*a.xx()*b.xx();
	c.d[1] = a.xx()*b.yy() + b.xx()*a.yy();
	c.d[3] = a.xx()*b.zz() + b.xx()*a.zz();
	c.d[6] = a.xx()*b.xy() + b.xx()*a.xy();
	c.d[10] = a.xx()*b.yz() + b.xx()*a.yz();
	c.d[15] = a.xx()*b.xz() + b.xx()*a.xz();

	c.d[2] = 2*a.yy()*b.yy();
	c.d[4] = a.yy()*b.zz() + b.yy()*a.zz();
	c.d[7] = a.yy()*b.xy() + b.yy()*a.xy();
	c.d[11] = a.yy()*b.yz() + b.yy()*a.yz();
	c.d[16] = a.yy()*b.xz() + b.yy()*a.xz();

	c.d[5] = 2*a.zz()*b.zz();
	c.d[8] = a.zz()*b.xy() + b.zz()*a.xy();
	c.d[12] = a.zz()*b.yz() + b.zz()*a.yz();
	c.d[17] = a.zz()*b.xz() + b.zz()*a.xz();

	c.d[9] = 2*a.xy()*b.xy();
	c.d[13] = a.xy()*b.yz() + b.xy()*a.yz();
	c.d[18] = a.xy()*b.xz() + b.xy()*a.xz();

	c.d[14] = 2*a.yz()*b.yz();
	c.d[19] = a.yz()*b.xz() + b.yz()*a.xz();

	c.d[20] = 2*a.xz()*b.xz();
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad4s a)_ijkl = (a_ik a_jl + a_il a_jk)/2
inline tens4ds dyad4s(const mat3ds& a)
{
	tens4ds c;
	c.d[0] =  a.xx()*a.xx(); // 0000 -> 00x00
	c.d[1] =  a.xy()*a.xy(); // 0011 -> 01x01
	c.d[3] =  a.xz()*a.xz(); // 0022 -> 02x02
	c.d[6] =  a.xx()*a.xy(); // 0001 -> 01x00
	c.d[10] = a.xy()*a.xz(); // 0012 -> 01x02
	c.d[15] = a.xx()*a.xz(); // 0002 -> 00x02

	c.d[2] =  a.yy()*a.yy(); // 1111 -> 11x11
	c.d[4] =  a.yz()*a.yz(); // 1122 -> 12x12
	c.d[7] =  a.xy()*a.yy(); // 1101 -> 11x01
	c.d[11] = a.yy()*a.yz(); // 1112 -> 11x12
	c.d[16] = a.xy()*a.yz(); // 1102 -> 12x01

	c.d[5] =  a.zz()*a.zz(); // 2222 -> 22x22
	c.d[8] =  a.xz()*a.yz(); // 2201 -> 12x02
	c.d[12] = a.yz()*a.zz(); // 2212 -> 22x12
	c.d[17] = a.xz()*a.zz(); // 2202 -> 22x02

	c.d[9] =  (a.xx()*a.yy() + a.xy()*a.xy())*0.5; // 0101 -> 0.5*(00x11 + 01x01)
	c.d[13] = (a.xy()*a.yz() + a.xz()*a.yy())*0.5; // 0112 -> 0.5*(01x12 + 02x11)
	c.d[18] = (a.xx()*a.yz() + a.xz()*a.xy())*0.5; // 0102 -> 0.5*(00x12 + 02x01)

	c.d[14] = (a.yy()*a.zz() + a.yz()*a.yz())*0.5; // 1212 -> 0.5*(11x22 + 12x12)
	c.d[19] = (a.xy()*a.zz() + a.xz()*a.yz())*0.5; // 1202 -> 0.5*(01x22 + 12x02)

	c.d[20] = (a.xx()*a.zz() + a.xz()*a.xz())*0.5; // 0202 -> 0.5*(00x22 + 02x02)

	return c;
}

//-----------------------------------------------------------------------------
// (a dyad4s b)_ijkl = (a_ik b_jl + a_il b_jk)/2 +  (b_ik a_jl + b_il a_jk)/2
inline tens4ds dyad4s(const mat3ds& a, const mat3ds& b)
{
	tens4ds c;
	c.d[0] =  2*a.xx()*b.xx();
	c.d[1] =  2*a.xy()*b.xy();
	c.d[3] =  2*a.xz()*b.xz();
	c.d[6] =  a.xy()*b.xx() + a.xx()*b.xy();
	c.d[10] = a.xz()*b.xy() + a.xy()*b.xz();
	c.d[15] = a.xz()*b.xx() + a.xx()*b.xz();

	c.d[2] =  2*a.yy()*b.yy();
	c.d[4] =  2*a.yz()*b.yz();
	c.d[7] =  a.yy()*b.xy() + a.xy()*b.yy();
	c.d[11] = a.yz()*b.yy() + a.yy()*b.yz();
	c.d[16] = a.yz()*b.xy() + a.xy()*b.yz();

	c.d[5] =  2*a.zz()*b.zz();
	c.d[8] =  a.yz()*b.xz() + a.xz()*b.yz();
	c.d[12] = a.zz()*b.yz() + a.yz()*b.zz();
	c.d[17] = a.zz()*b.xz() + a.xz()*b.zz();

	c.d[9] =  0.5*(a.yy()*b.xx() + 2*a.xy()*b.xy() + a.xx()*b.yy());
	c.d[13] = 0.5*(a.yz()*b.xy() + a.yy()*b.xz() + a.xz()*b.yy() + a.xy()*b.yz());
	c.d[18] = 0.5*(a.yz()*b.xx() + a.xz()*b.xy() + a.xy()*b.xz() + a.xx()*b.yz());

	c.d[14] = 0.5*(a.zz()*b.yy() + 2*a.yz()*b.yz() + a.yy()*b.zz());
	c.d[19] = 0.5*(a.zz()*b.xy() + a.yz()*b.xz() + a.xz()*b.yz() + a.xy()*b.zz());

	c.d[20] = 0.5*(a.zz()*b.xx() + 2*a.xz()*b.xz() + a.xx()*b.zz());
	return c;

}
