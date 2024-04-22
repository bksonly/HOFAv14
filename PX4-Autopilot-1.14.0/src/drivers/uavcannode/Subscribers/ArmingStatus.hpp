/****************************************************************************
 *
 *   Copyright (c) 2021 PX4 Development Team. All rights reserved.
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

#pragma once

#include "UavcanSubscriberBase.hpp"

#include <uavcan/equipment/safety/ArmingStatus.hpp>

#include <uORB/Publication.hpp>
#include <uORB/topics/actuator_armed.h>

namespace uavcannode
{

class ArmingStatus;

typedef uavcan::MethodBinder<ArmingStatus *,
	void (ArmingStatus::*)(const uavcan::ReceivedDataStructure<uavcan::equipment::safety::ArmingStatus>&)>
	ArmingStatusBinder;

class ArmingStatus :
	public UavcanSubscriberBase,
	private uavcan::Subscriber<uavcan::equipment::safety::ArmingStatus, ArmingStatusBinder>
{
public:
	ArmingStatus(uavcan::INode &node) :
		UavcanSubscriberBase(uavcan::equipment::safety::ArmingStatus::DefaultDataTypeID),
		uavcan::Subscriber<uavcan::equipment::safety::ArmingStatus, ArmingStatusBinder>(node)
	{}

	bool init()
	{
		if (start(ArmingStatusBinder(this, &ArmingStatus::callback)) < 0) {
			PX4_ERR("uavcan::equipment::safety::ArmingStatus subscription failed");
			return false;
		}

		return true;
	}

	void PrintInfo() const override
	{
		printf("\t%s:%d -> %s\n",
		       uavcan::equipment::safety::ArmingStatus::getDataTypeFullName(),
		       uavcan::equipment::safety::ArmingStatus::DefaultDataTypeID,
		       _actuator_armed_pub.get_topic()->o_name);
	}

private:
	void callback(const uavcan::ReceivedDataStructure<uavcan::equipment::safety::ArmingStatus> &msg)
	{
		actuator_armed_s actuator_armed{};

		if (msg.status) {
			actuator_armed.armed = true;

		} else {
			actuator_armed.armed = false;
		}

		actuator_armed.timestamp = hrt_absolute_time();
		_actuator_armed_pub.publish(actuator_armed);
	}

	uORB::Publication<actuator_armed_s> _actuator_armed_pub{ORB_ID(actuator_armed)};
};
} // namespace uavcannode
