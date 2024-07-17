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
 * @file RateControl.cpp
 */

#include "rate_control.hpp"
#include <px4_platform_common/defines.h>

using namespace matrix;

void RateControl::setGains(const Vector3f &P, const Vector3f &I, const Vector3f &D)
{
	_gain_p = P;
	_gain_i = I;
	_gain_d = D;
}

void RateControl::setSaturationStatus(const Vector3<bool> &saturation_positive,
				      const Vector3<bool> &saturation_negative)
{
	_control_allocator_saturation_positive = saturation_positive;
	_control_allocator_saturation_negative = saturation_negative;
}

void RateControl::setPositiveSaturationFlag(size_t axis, bool is_saturated)
{
	if (axis < 3) {
		_control_allocator_saturation_positive(axis) = is_saturated;
	}
}

void RateControl::setNegativeSaturationFlag(size_t axis, bool is_saturated)
{
	if (axis < 3) {
		_control_allocator_saturation_negative(axis) = is_saturated;
	}
}

Vector3f RateControl::update(const Vector3f &rate, const Vector3f &rate_sp, const Vector3f &angular_accel,
			     const float dt, const bool landed)
{
	// angular rates error
	Vector3f rate_error = rate_sp - rate;

	// PID control with feed forward
	const Vector3f torque = _gain_p.emult(rate_error) + _rate_int - _gain_d.emult(angular_accel) + _gain_ff.emult(rate_sp);

	// update integral only if we are not landed
	if (!landed) {
		updateIntegral(rate_error, dt);
	}

	return torque;
}

Vector3f RateControl::update(const Quatf &q,const Vector3f &rate, const Vector3f &rate_sp, const Vector3f &angular_accel,const float dt, const bool landed)
{
	Dcmf Rd(_attitude_setpoint_q);
	Dcmf R(q);
	Dcmf RdT = Rd.transpose();
	Dcmf RT = R.transpose();

	Vector3f w = rate;
	Vector3f wd = rate_sp;
	Dcmf W=w.hat();
	Dcmf Wd=wd.hat();

	Vector3f e;
	Vector3f de;
	Dcmf ehat;

	ehat = RdT * R - RT * Rd;
	e = ehat.vee();
	de = (e - pre_e)/dt;
	pre_e = e;

	Vector3f ew = w - RT*Rd*wd;

	Dcmf tmp = RdT*R*ew.hat()*W - Wd*RdT*R*ew.hat();
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
	J(0,0) = 0.029125;
	J(1,1) = 0.029125;
	J(2,2) = 0.055225;

	SquareMatrix<float, 3> B;
	SquareMatrix<float, 3> invB;
	B = RB * J.I();
	invB = B.I();

	Vector3f M_star;
	SquareMatrix<float, 3> Ke;
	SquareMatrix<float, 3> Kde;
	Ke(0,0) = 10;
	Ke(1,1) = 10;
	Ke(2,2) = 7.07;
	Kde(0,0) = 4.79;
	Kde(1,1) = 4.79;
	Kde(2,2) = 3.89;
	M_star = -Ke * e -Kde * de;

	Vector3f M;

	
	M = -invB * Ad + w.cross(J*w) +invB * M_star;

	M(0) = math::constrain(M(0), -0.5f, 0.5f);
	M(1) = math::constrain(M(1), -0.5f, 0.5f);
	M(2) = math::constrain(M(2), -0.04f, 0.04f);
	return M;
	// // angular rates error
	// Vector3f rate_error = rate_sp - rate;

	// // PID control with feed forward
	// const Vector3f torque = _gain_p.emult(rate_error) + _rate_int - _gain_d.emult(angular_accel) + _gain_ff.emult(rate_sp);

	// // update integral only if we are not landed
	// if (!landed) {
	// 	updateIntegral(rate_error, dt);
	// }

	// return torque;
}

void RateControl::updateIntegral(Vector3f &rate_error, const float dt)
{
	for (int i = 0; i < 3; i++) {
		// prevent further positive control saturation
		if (_control_allocator_saturation_positive(i)) {
			rate_error(i) = math::min(rate_error(i), 0.f);
		}

		// prevent further negative control saturation
		if (_control_allocator_saturation_negative(i)) {
			rate_error(i) = math::max(rate_error(i), 0.f);
		}

		// I term factor: reduce the I gain with increasing rate error.
		// This counteracts a non-linear effect where the integral builds up quickly upon a large setpoint
		// change (noticeable in a bounce-back effect after a flip).
		// The formula leads to a gradual decrease w/o steps, while only affecting the cases where it should:
		// with the parameter set to 400 degrees, up to 100 deg rate error, i_factor is almost 1 (having no effect),
		// and up to 200 deg error leads to <25% reduction of I.
		float i_factor = rate_error(i) / math::radians(400.f);
		i_factor = math::max(0.0f, 1.f - i_factor * i_factor);

		// Perform the integration using a first order method
		float rate_i = _rate_int(i) + i_factor * _gain_i(i) * rate_error(i) * dt;

		// do not propagate the result if out of range or invalid
		if (PX4_ISFINITE(rate_i)) {
			_rate_int(i) = math::constrain(rate_i, -_lim_int(i), _lim_int(i));
		}
	}
}

void RateControl::getRateControlStatus(rate_ctrl_status_s &rate_ctrl_status)
{
	rate_ctrl_status.rollspeed_integ = _rate_int(0);
	rate_ctrl_status.pitchspeed_integ = _rate_int(1);
	rate_ctrl_status.yawspeed_integ = _rate_int(2);
}
