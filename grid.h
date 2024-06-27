#pragma once

template<int GRID_SIZE>
struct grid_t {
	static const bool use_blocked_layout = false;

	float values[GRID_SIZE * GRID_SIZE];

	grid_t(float init_v) {
        for (auto& v : values) {
            v = init_v;
        }
	}

	float& at_ref(int x, int y) {
		if (use_blocked_layout) {
			static const int BLOCK_SIZE = 4;
			static_assert(GRID_SIZE % BLOCK_SIZE == 0, "layout in 4x4 blocks");

			int b_x = x / BLOCK_SIZE;
			int b_y = y / BLOCK_SIZE;
			// int in_b_x = x - b_x * BLOCK_SIZE;
			// int in_b_y = y - b_y * BLOCK_SIZE;
			// auto block_offset = (b_y * GRID_SIZE + b_x * BLOCK_SIZE) * BLOCK_SIZE;
			// return values.data()[block_offset + in_b_y * BLOCK_SIZE + in_b_x];

			auto data_offset = (b_y * (GRID_SIZE - BLOCK_SIZE) + b_x * (BLOCK_SIZE - 1) + y) * BLOCK_SIZE + x;
			return values[data_offset];
		} else {
			return values[y * GRID_SIZE + x];
		}
	}

	float at(int x, int y) {
		return at_ref(x < GRID_SIZE ? x : GRID_SIZE - 1, y < GRID_SIZE ? y : GRID_SIZE - 1);
	}

	void set(int x, int y, float v) {
		at_ref(x, y) = v;
	}
};
