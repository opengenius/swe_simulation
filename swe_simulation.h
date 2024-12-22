#pragma once

#include "grid.h"

// Real-time Simulation of Large Bodies of Water with Small Scale Details
// https://matthias-research.github.io/pages/publications/hfFluid.pdf

template<int GRID_SIZE>
struct swe_simulation_t {
    using grid_sized_t = grid_t<GRID_SIZE>;

    const float dxdy = 0.4f;

    grid_sized_t water_heights[2];
	grid_sized_t vel_us[2];
	grid_sized_t vel_vs[2];
	grid_sized_t ground_h;
	int water_index = 0;

	void set_height(int x, int y, float v) {
		water_heights[water_index].set(x, y, v);
	}

	float height(int x, int y) const {
		return water_heights[water_index].at(x, y);
	}

	float velocity_u(int x, int y) const {
		return vel_us[water_index].at(x, y, 0.0f);
	}

	float velocity_v(int x, int y) const {
		return vel_vs[water_index].at(x, y, 0.0f);
	}

	static float interpolate(float a, float b, float w) {
		return a * (1.0f - w) + b * w;
	}

	static float interpolate_grid(const grid_sized_t& grid, float x, float y) {
		int xi = (int)x;
		int yi = (int)y;
		float x_fr = x - xi;
		float y_fr = y - yi;

		float g00 = grid.at(xi, yi, 0.0f);
		float g10 = grid.at(xi + 1, yi, 0.0f);
		float g01 = grid.at(xi, yi + 1, 0.0f);
		float g11 = grid.at(xi + 1, yi + 1, 0.0f);

		return interpolate(
			interpolate(g00, g10, x_fr),
			interpolate(g01, g11, x_fr),
			y_fr);
	}

    void step() {
        const float g = 9.81f;
		const float dt = 0.01666f;	
		const float EPS = 0.001f;//0.0001f * dxdy;
		const float max_vel = dxdy / dt * 0.5f; // 0.25f seems more stable with lower dx
        const float drag_shore_height_threshold = 0.1f;
        const float drag_factor = 0.0f;

		int next_water_index = (water_index + 1) % 2;
		
		// advect velocities (Semi-Lagrangian Method)
		auto& vel_u_src = vel_us[water_index];
		auto& vel_v_src = vel_vs[water_index];
		auto& vel_u_dst = vel_us[next_water_index];
		auto& vel_v_dst = vel_vs[next_water_index];
		for (int yi = 0; yi < GRID_SIZE - 1; ++yi) {
			for (int xi = 0; xi < GRID_SIZE - 1; ++xi) {
				{
					auto u_xy = vel_u_src.at(xi, yi);
					auto v_xy = interpolate_grid(vel_v_src, xi + 0.5f, yi - 0.5f);

					float pos[2] = {float(xi), float(yi)};
					pos[0] -= dt * u_xy;
					pos[1] -= dt * v_xy;

					float new_u = interpolate_grid(vel_u_src, pos[0], pos[1]);

					vel_u_dst.set(xi, yi, new_u);
				}

				{
					auto v_xy = vel_v_src.at(xi, yi);
					auto u_xy = interpolate_grid(vel_u_src, xi - 0.5f, yi + 0.5f);

					float pos[2] = {float(xi), float(yi)};
					pos[0] -= dt * u_xy;
					pos[1] -= dt * v_xy;

					float new_v = interpolate_grid(vel_v_src, pos[0], pos[1]);

					vel_v_dst.set(xi, yi, new_v);
				}
				
			}
		}

        // update velocities
		auto& water_h = water_heights[water_index];
		auto& water_h_dst = water_heights[next_water_index];
		auto& vel_u = vel_us[next_water_index];
		auto& vel_v = vel_vs[next_water_index];
		for (int yi = 0; yi < GRID_SIZE - 1; ++yi) {
			for (int xi = 0; xi < GRID_SIZE - 1; ++xi) {
				float h_ij = water_h.at_clamped(xi, yi);
				float h_i1j = water_h.at_clamped(xi + 1, yi);
				float h_ij1 = water_h.at_clamped(xi, yi + 1);

				float H_ij = ground_h.at_clamped(xi, yi);
				float H_i1j = ground_h.at_clamped(xi + 1, yi);
				float H_ij1 = ground_h.at_clamped(xi, yi + 1);

				float n_ij = h_ij + H_ij;
				float n_i1j = h_i1j + H_i1j;
				float n_ij1 = h_ij1 + H_ij1;

				float u_ref = vel_u.at(xi, yi);
				float v_ref = vel_v.at(xi, yi);

				// enable drag force on the shores
				float drag_koef = h_ij < drag_shore_height_threshold ? drag_factor : 0.0f;

				if ((n_i1j < H_ij && h_ij < EPS) ||
					n_ij < H_i1j && h_i1j < EPS) {
					u_ref = 0.0f;
				} else {
					float dh_dx = (n_i1j - n_ij) / dxdy;
					float new_u = u_ref - dt * g * dh_dx - drag_koef * u_ref;
					u_ref = std::max(-max_vel, std::min(new_u, max_vel));
				}

				if ((n_ij1 < H_ij && h_ij < EPS) ||
					n_ij < H_ij1 && h_ij1 < EPS) {
					v_ref = 0.0f;
				} else {
					float dh_dy = (n_ij1 - n_ij) / dxdy;
					float new_v = v_ref - dt * g * dh_dy - drag_koef * v_ref;
					v_ref = std::max(-max_vel, std::min(new_v, max_vel));
				}

				vel_u.set(xi, yi, u_ref);
				vel_v.set(xi, yi, v_ref);
			}
		}
		static const bool limit_heights = false;
		
		// update height
		for (int yi = 0; yi < GRID_SIZE; ++yi) {
			for (int xi = 0; xi < GRID_SIZE; ++xi) {
				float u_ij = vel_u.at(xi, yi);
				float u_im1j = 0 < xi ? vel_u.at(xi - 1, yi) : 0.0f;
				float v_ij = vel_v.at(xi, yi);
				float v_ijm1 = 0 < yi ? vel_v.at(xi, yi - 1) : 0.0f;
				float h_ij = water_h.at_clamped(xi, yi);

				float h_ip2_j = u_ij   <= 0.0f ? water_h.at_clamped(xi + 1, yi) : h_ij;
				float h_im2_j = u_im1j <= 0.0f ? h_ij : water_h.at_clamped(xi - 1, yi);
				float h_i_jp2 = v_ij   <= 0.0f ? water_h.at_clamped(xi, yi + 1) : h_ij;
				float h_i_jm2 = v_ijm1 <= 0.0f ? h_ij : water_h.at_clamped(xi, yi - 1);

				if (limit_heights) {
					const float b = 2.0f;
					float h_avgmax = b * dxdy / (g * dt);
					float h_adj = std::max(0.0f, (h_ip2_j + h_im2_j + h_i_jp2 + h_i_jm2) * 0.25f - h_avgmax);
					h_ip2_j -= h_adj;
					h_im2_j -= h_adj;
					h_i_jp2 -= h_adj;
					h_i_jm2 -= h_adj;
				}

				float dh_dt = (h_ip2_j * u_ij - h_im2_j * u_im1j) / dxdy +
							(h_i_jp2 * v_ij - h_i_jm2 * v_ijm1) / dxdy;

				float new_h = std::max(0.0f, h_ij - dh_dt * dt);
				assert(std::isfinite(new_h));

				water_h_dst.set(xi, yi, new_h);
			}
		}

		water_index = next_water_index;
    }

    float get_divergence(int xi, int yi) {
		float u11 = velocity_u(xi, yi);
		float u01 = velocity_u(xi - 1, yi);
		float v11 = velocity_v(xi, yi);
		float v10 = velocity_v(xi, yi - 1);

        return (u11 - u01 + v11 - v10) / dxdy;
    }
};