/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file AttitudeControl.cpp
 */

#include <AttitudeControl.hpp>
#include <mathlib/math/Functions.hpp>

using namespace matrix;


void AttitudeControl::setProportionalGain(const matrix::Vector3f &proportional_gain, const float yaw_weight)
{
	_proportional_gain = proportional_gain;
	_yaw_w = math::constrain(yaw_weight, 0.f, 1.f);

	// compensate for the effect of the yaw weight rescaling the output
	if (_yaw_w > 1e-4f) {
		_proportional_gain(2) /= _yaw_w;
	}
}

matrix::Vector3f AttitudeControl::update(const Quatf &q, const matrix::Vector3f &rates, const matrix::Vector3f &angular_accel, const float dt_R)
{

	
	const float dt_setRd = math::constrain(((_current_setRd_timestamp - _previous_setRd_timestamp) * 1e-6f), 0.000125f, 0.02f);//其实根本不知道上层更新频率，大概就这样罢，这个dt实质上也不是两个R的更新间隔，只是这一层的读取间隔，干脆也作为角速度的更新间隔 改后 这下精确了 和Rd同时产生的

	//_previous_updateR_timestamp = _current_updateR_timestamp;
	//_current_updateR_timestamp = hrt_absolute_time(); // Update the current timestamp 1000000us=1s

	Quatf qd = _attitude_setpoint_q;
	Quatf pre_qd = _previous_attitude_setpoint_q;

	//const float dt_updateR = math::constrain(((_current_updateR_timestamp -       _previous_updateR_timestamp) * 1e-6f), 0.000125f, 0.02f);

	const float dt_updateR = dt_R;//精确的，限幅在mc_att_control_main里

	Dcmf Rd(qd);
	Dcmf R(q);
	Dcmf RdT = Rd.transpose();
	Dcmf RT = R.transpose();

	Vector3f e;
	Vector3f de;
	Dcmf ehat;
	ehat = RdT * R - RT * Rd;
	e = ehat.vee();
	de = (e - _pre_e)/dt_updateR;
	_pre_e = e;

	if (pre_qd.dot(qd) < 0) {
    qd = -qd; // 取反四元数，保证正负一致性
	}


	Quatf Q_wd = qd.inversed() * (qd - pre_qd) / dt_setRd * 2.0f;
	Vector3f rate_setpoint = Q_wd.imag();

	if (_pre_q.dot(q) < 0) {
    _pre_q = -_pre_q; // 取反四元数，保证正负一致性
	}
	

	//Quatf Q_w = q.inversed() * (q - _pre_q) / dt_updateR * 2.0f;
	//Vector3f w = Q_w.imag();

	Vector3f w = rates;

	for (int i = 0; i < 3; i++) {
		rate_setpoint(i) = math::constrain(rate_setpoint(i), -_rate_limit(i), _rate_limit(i));
	}

	for (int i = 0; i < 3; i++) {
		w(i) = math::constrain(w(i), -_rate_limit(i), _rate_limit(i));
	}


	Vector3f wd = rate_setpoint;
	Vector3f ad = (wd - _pre_rate_setpoint)/dt_setRd;
	_pre_rate_setpoint = rate_setpoint;

	for (int i = 0; i < 3; i++) {
		ad(i) = math::constrain(ad(i), -10.f, 10.f);
	}

	Dcmf W=w.hat();
	Dcmf Wd=wd.hat();
	Vector3f ew = w - RT*Rd*wd;

	Dcmf tmp = RdT*R*ew.hat()*W - Wd*RdT*R*ew.hat() - ad.hat()*RdT*R;
	Dcmf Adhat = (tmp - tmp.transpose());
	Vector3f Ad = Adhat.vee();

	Dcmf RB;
	Dcmf Rr = RdT * R;

	RB(0, 0) = Rr(1, 1) + Rr(2, 2);
	RB(0, 1) = -Rr(1, 0);
	RB(0, 2) = -Rr(2, 0);

	RB(1, 0) = -Rr(0, 1);
	RB(1, 1) = Rr(0, 0) + Rr(2, 2);
	RB(1, 2) = -Rr(2, 1);

	RB(2, 0) = -Rr(0, 2);
	RB(2, 1) = -Rr(1, 2);
	RB(2, 2) = Rr(0, 0) + Rr(1, 1);
	
	SquareMatrix<float, 3> J;
	J.setZero();
	J(0,0) = 0.0024;
	J(1,1) = 0.0025;
	J(2,2) = 0.0041;

	SquareMatrix<float, 3> B;
	SquareMatrix<float, 3> invB;
	B = RB * J.I();

	Vector3f M_star;
	SquareMatrix<float, 3> Ke;
	SquareMatrix<float, 3> Kde;
	Ke(0,0) = 10;
	Ke(1,1) = 10;
	Ke(2,2) = 3.1622;
	Kde(0,0) = 5.477;
	Kde(1,1) = 5.477;
	Kde(2,2) = 4.040;
	M_star = -Ke * e -Kde * de;

	Vector3f M;

	bool is_invertible = matrix::inv(B, invB);

    if (is_invertible) {
        	M = -invB * Ad + invB * M_star;
    } else {
		M(0)=0;
		M(1)=0;
		M(2)=0;
    }
	//return rate_setpoint;
	return M;//只能通过这个通道了，wxJw还没加上
}
