from typing import TYPE_CHECKING

from vulkan import (
    VK_API_VERSION_1_0,
    VK_MAKE_VERSION,
    VK_STRUCTURE_TYPE_APPLICATION_INFO,
    VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
    VK_VERSION_MAJOR,
    VK_VERSION_MINOR,
    VK_VERSION_PATCH,
    VkApplicationInfo,
    VkError,
    VkInstanceCreateInfo,
    vkCreateInstance,
    vkDestroyInstance,
    vkEnumeratePhysicalDevices,
    vkGetPhysicalDeviceProperties,
)

if TYPE_CHECKING:
    from typing import Optional

DEVICE_ID = "nvidia"


def get_available_gpu_devices():
    app_info = VkApplicationInfo(
        sType=VK_STRUCTURE_TYPE_APPLICATION_INFO,
        pApplicationName=b"VulkanChecker",
        applicationVersion=VK_MAKE_VERSION(1, 0, 0),
        pEngineName=b"NoEngine",
        engineVersion=VK_MAKE_VERSION(1, 0, 0),
        apiVersion=VK_API_VERSION_1_0,
    )

    instance_info = VkInstanceCreateInfo(
        sType=VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO, pApplicationInfo=app_info
    )

    try:
        instance = vkCreateInstance(instance_info, None)
    except VkError as e:
        print(f"Failed to create Vulkan instance: {e}")
        return

    devices = vkEnumeratePhysicalDevices(instance)
    devices_data = []

    for i, device in enumerate(devices):
        props = vkGetPhysicalDeviceProperties(device)
        device_name = props.deviceName
        device_id = i
        device_type = {
            1: "1 (Integrated GPU)",
            2: "2 (Discrete GPU)",
        }.get(props.deviceType, f"{props.deviceType} (Other)")
        api_version = props.apiVersion

        devices_data.append(
            {
                "device_name": device_name,
                "device_id": device_id,
                "device_type": device_type,
                "vulkan_version": f"{VK_VERSION_MAJOR(api_version)}.{VK_VERSION_MINOR(api_version)}.{VK_VERSION_PATCH(api_version)}",
            }
        )

    vkDestroyInstance(instance, None)

    return devices_data


def get_selected_gpu_device(
    devices_data: "Optional[list[dict[str, str]]]", device_id: str | int = DEVICE_ID
) -> "Optional[dict[str, str]]":
    if not devices_data:
        devices_data = get_available_gpu_devices()

    for device in devices_data:
        name = device["device_name"].strip().lower()
        i = device["device_id"]

        if isinstance(device_id, str) and device_id in name:
            return device

        elif isinstance(device_id, int) and device_id == i:
            return device

    return devices_data[0] if devices_data else None


def list_available_gpu_devices():
    devices_data = get_available_gpu_devices()

    if not devices_data:
        print("No Vulkan devices found.")
        return

    selected_device = get_selected_gpu_device(devices_data)

    for device in devices_data:
        name = device["device_name"]
        device_id = device["device_id"]
        device_type = device["device_type"]
        vulkan_version = device["vulkan_version"]

        s = "✔" if device == selected_device else " "

        print(f"[{s}] Device Name    : {name}")
        print(f"    Device ID      : {device_id}")
        print(f"    Device Type    : {device_type}")
        print(f"    Vulkan Version : {vulkan_version}")
        print()
