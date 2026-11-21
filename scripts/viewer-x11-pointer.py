#!/usr/bin/env python3
"""Deliver the native X11 drag and scroll for the Sirius viewer check."""

import ctypes as c
import sys
import time


def deliver(window):
    x11 = c.CDLL("libX11.so.6")
    xtest = c.CDLL("libXtst.so.6")
    display_type, window_type = c.c_void_p, c.c_ulong
    integer, pointer = c.c_int, c.POINTER(c.c_int)
    window_pointer = c.POINTER(window_type)
    x11.XOpenDisplay.argtypes = (c.c_char_p,)
    x11.XOpenDisplay.restype = display_type
    x11.XDefaultRootWindow.argtypes = (display_type,)
    x11.XDefaultRootWindow.restype = window_type
    x11.XDefaultScreen.argtypes = (display_type,)
    for dimension in (x11.XDisplayWidth, x11.XDisplayHeight):
        dimension.argtypes = (display_type, integer)
    x11.XTranslateCoordinates.argtypes = (
        display_type, window_type, window_type, integer, integer, pointer, pointer, window_pointer,
    )
    x11.XQueryPointer.argtypes = (
        display_type, window_type, window_pointer, window_pointer,
        pointer, pointer, pointer, pointer, c.POINTER(c.c_uint),
    )
    x11.XSync.argtypes = (display_type, integer)
    x11.XCloseDisplay.argtypes = (display_type,)
    xtest.XTestQueryExtension.argtypes = (display_type, pointer, pointer, pointer, pointer)
    xtest.XTestFakeMotionEvent.argtypes = (display_type, integer, integer, integer, c.c_ulong)
    xtest.XTestFakeButtonEvent.argtypes = (display_type, c.c_uint, integer, c.c_ulong)

    display = x11.XOpenDisplay(None)
    if not display:
        raise RuntimeError("could not open the viewer's X11 display")
    try:
        if not xtest.XTestQueryExtension(display, *(c.byref(integer()) for _ in range(4))):
            raise RuntimeError("the viewer's X11 display lacks XTEST")
        root = x11.XDefaultRootWindow(display)
        screen = x11.XDefaultScreen(display)
        width, height = (x11.XDisplayWidth(display, screen), x11.XDisplayHeight(display, screen))

        def translate(x, y):
            sx, sy, child = integer(), integer(), window_type()
            if not x11.XTranslateCoordinates(
                display, window, root, x, y, c.byref(sx), c.byref(sy), c.byref(child),
            ) or not (0 <= sx.value < width and 0 <= sy.value < height):
                raise RuntimeError("the viewer input target is outside the visible X11 screen")
            return sx.value, sy.value

        start, end = translate(120, 120), translate(180, 160)

        def move(target):
            if not xtest.XTestFakeMotionEvent(display, screen, *target, 0):
                raise RuntimeError("XTEST could not move the pointer")
            x11.XSync(display, 0)
            # Walk the visible window tree at the actual pointer. A covered
            # or unmapped viewer must never send clicks to another window.
            child = root
            while child and child != window:
                actual_root, below = window_type(), window_type()
                coordinates = [integer() for _ in range(4)]
                mask = c.c_uint()
                if not x11.XQueryPointer(
                    display, child, c.byref(actual_root), c.byref(below),
                    *(c.byref(value) for value in coordinates), c.byref(mask),
                ):
                    raise RuntimeError("could not locate the pointer after XTEST motion")
                child = below.value
            if child != window:
                raise RuntimeError("the viewer is covered or has no visible input surface")
            time.sleep(0.1)

        def button(index, pressed):
            if not xtest.XTestFakeButtonEvent(display, index, pressed, 0):
                raise RuntimeError("XTEST could not deliver a pointer button")
            x11.XSync(display, 0)
            time.sleep(0.1)

        move(start)
        button(1, True)
        try:
            move(end)
        finally:
            button(1, False)
        button(4, True)
        button(4, False)
    finally:
        x11.XCloseDisplay(display)


if __name__ == "__main__":
    try:
        deliver(int(sys.argv[1]))
    except (IndexError, ValueError, OSError, RuntimeError) as error:
        raise SystemExit(f"Native viewer pointer input failed: {error}") from None
