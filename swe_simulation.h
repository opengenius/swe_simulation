#pragma once

#include "grid.h"

// Real-time Simulation of Large Bodies of Water with Small Scale Details
// https://matthias-research.github.io/pages/publications/hfFluid.pdf

template<int GRID_SIZE>
struct swe_simulation_t {
    using grid_sized_t = grid_t<GRID_SIZE>;

    grid_sized_t water_h = grid_sized_t(0.0f);
	grid_sized_t water_h_dst = grid_sized_t(0.0f);
	grid_sized_t vel_u = grid_sized_t(0.0f);
	grid_sized_t vel_v = grid_sized_t(0.0f);
	grid_sized_t ground_h = grid_sized_t(0.0f);

    const float dxdy = 0.4f;
    void step() {
        const float g = 9.81f;
		const float dt = 0.01666f;	
		const float EPS = 0.0001f * dxdy;
		const float max_vel = dxdy / dt * 0.5f;
        const float drag_shore_height_threshold = 0.1f;
        const float drag_factor = 0.05f;

        // update velocities
		for (int yi = 0; yi < GRID_SIZE - 1; ++yi) {
			for (int xi = 0; xi < GRID_SIZE - 1; ++xi) {
				float h_ij = water_h.at(xi, yi);
				float h_i1j = water_h.at(xi + 1, yi);
				float h_ij1 = water_h.at(xi, yi + 1);

				float H_ij = ground_h.at(xi, yi);
				float H_i1j = ground_h.at(xi + 1, yi);
				float H_ij1 = ground_h.at(xi, yi + 1);

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

		// update height
		for (int yi = 0; yi < GRID_SIZE; ++yi) {
			for (int xi = 0; xi < GRID_SIZE; ++xi) {
				float u_ij = vel_u.at(xi, yi);
				float u_im1j = 0 < xi ? vel_u.at(xi - 1, yi) : 0.0f;
				float v_ij = vel_v.at(xi, yi);
				float v_ijm1 = 0 < yi ? vel_v.at(xi, yi - 1) : 0.0f;
				float h_ij = water_h.at(xi, yi);

				float h_ip2_j = u_ij   <= 0.0f ? water_h.at(xi + 1, yi) : h_ij;
				float h_im2_j = u_im1j <= 0.0f ? h_ij     			    : water_h.at(xi - 1, yi);
				float h_i_jp2 = v_ij   <= 0.0f ? water_h.at(xi, yi + 1) : h_ij;
				float h_i_jm2 = v_ijm1 <= 0.0f ? h_ij     			    : water_h.at(xi, yi - 1);

				float dh_dt = (h_ip2_j * u_ij - h_im2_j * u_im1j) / dxdy +
							(h_i_jp2 * v_ij - h_i_jm2 * v_ijm1) / dxdy;

				float new_h = std::max(0.0f, h_ij - dh_dt * dt);
				assert(std::isfinite(new_h));

				water_h_dst.set(xi, yi, new_h);
			}
		}
        std::swap(water_h, water_h_dst);
    }

    float get_divergence(int xi, int yi) {
        // todo: is it a divergence?)

        return (vel_u.at(xi, yi) - vel_u.at(xi - 1, yi)) / dxdy + 
            (vel_v.at(xi, yi) - vel_v.at(xi, yi - 1)) / dxdy;
    }
};