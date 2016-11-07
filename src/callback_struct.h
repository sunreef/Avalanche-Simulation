#pragma once

enum ModifierKeys {
	SHIFT, CTRL, ALT
};

struct CallbackStruct {
	int x, y;
	int dx, dy;

	bool ascii_keys[256];
	bool mouse_buttons[3];
	bool modifier_keys[3];

	CallbackStruct() {
		x = 0;
		y = 0;
		dx = 0;
		dy = 0;

		for (int i = 0; i < 256; i++) {
			ascii_keys[i] = false;
		}
		for (int i = 0; i < 3; i++) {
			mouse_buttons[i] = false;
		}
		modifier_keys[0] = false;
		modifier_keys[1] = false;
		modifier_keys[2] = false;
	}
};
