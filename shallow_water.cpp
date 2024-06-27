#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

// Include glfw3.h after our OpenGL definitions
#include <GLFW/glfw3.h>
#include <cmath>
#include <chrono>

#include "grid.h"
#include "swe_simulation.h"

int main() {
	const int GRID_SIZE = 80;
	const float cell_size = 16.0f;

	int w = int(cell_size * GRID_SIZE);
	int h = int(cell_size * GRID_SIZE);

    if (!glfwInit())
        return -1;

    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);      

    GLFWwindow *window = glfwCreateWindow(w, h, "Shallow Water Simulation", NULL, NULL);
    if (window == NULL)
        return -1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

	swe_simulation_t<GRID_SIZE> sim = {};

	static const bool init_terrain_with_deepening = true;

	if (init_terrain_with_deepening) {
		// init test height map
		float center_x = GRID_SIZE / 2;
		float center_y = GRID_SIZE / 2;
		float radius = GRID_SIZE * 0.5f - 1.0f;
		float radius2 = radius * radius;
		for (int yi = 0; yi < GRID_SIZE; ++yi) {
			for (int xi = 0; xi < GRID_SIZE; ++xi) {
				float dx = center_x - (xi + 0.5f); // 0.5 is cell center
				float dy = center_y - (yi + 0.5f); // 0.5 is cell center
				float dist2 = dx * dx + dy * dy;

				sim.ground_h.set(xi, yi, sqrt(dist2) / radius);
				auto h = std::max(0.0f, 1.0f - sim.ground_h.at(xi, yi));
				// water_h_dst.set(xi, yi, h);
				// water_h.set(xi, yi, h);
			}
		}
	}

	const float cell_11_size = float(GRID_SIZE) / 11;
	const float cell_5_size = float(GRID_SIZE) / 5;
	for (int yi = 0; yi < GRID_SIZE; ++yi) {
		for (int xi = 0; xi < GRID_SIZE; ++xi) {
			if (xi < GRID_SIZE / 2 - 1) {
				float wh = std::max(0.0f, 1.0f - sim.ground_h.at(xi, yi));
				sim.water_h_dst.set(xi, yi, wh);
				sim.water_h.set(xi, yi, wh);
			}

			// 22 / 11 = 2
			if (xi == GRID_SIZE / 2) {
				auto cell11_index = int(float(yi) / cell_11_size);
				if ((cell11_index + 1) % 3 == 0) {
					// keep 0
				} else {
					sim.ground_h.set(xi, yi, 2.0f);
				}
			}

			if (xi == (GRID_SIZE * 3 / 4)) {
				auto cell5_index = int(float(yi) / cell_5_size);
				if (cell5_index == 2) {
					// keep 0
				} else {
					sim.ground_h.set(xi, yi, 2.0f);
				}
			}
		}
	}
	
	while (!glfwWindowShouldClose(window))
	{
		glfwPollEvents();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();

		ImGui::NewFrame();

		const ImGuiViewport* viewport = ImGui::GetMainViewport();
		const bool use_work_area = true;
		ImGui::SetNextWindowPos(use_work_area ? viewport->WorkPos : viewport->Pos);
		ImGui::SetNextWindowSize(use_work_area ? viewport->WorkSize : viewport->Size);

		ImGuiWindowFlags window_flags = 0;
		window_flags |= ImGuiWindowFlags_NoTitleBar;
		window_flags |= ImGuiWindowFlags_NoMove;
		window_flags |= ImGuiWindowFlags_NoResize;
		window_flags |= ImGuiWindowFlags_NoCollapse;
		window_flags |= ImGuiWindowFlags_NoBackground;
		window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus;

		ImGui::Begin("Simulation", nullptr, window_flags);
		
		// mouse interaction
		if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
			auto pos = ImGui::GetMousePos();
			pos.x = std::max(0.0f, pos.x / cell_size);
			pos.y = std::max(0.0f, pos.y / cell_size);

			int x = int(pos.x) < GRID_SIZE - 1 ? int(pos.x) : GRID_SIZE - 1;
			int y = int(pos.y) < GRID_SIZE - 1 ? int(pos.y) : GRID_SIZE - 1;

			if (x != 0 && y != 0 && x != (GRID_SIZE - 1) && y != (GRID_SIZE - 1)) {
				sim.water_h.set(x, y, 1.1f);
				// vel_u.at(x, y) = 10.0f;
			}
		}

		auto start = std::chrono::high_resolution_clock::now();
		sim.step();

		auto elapsed = std::chrono::high_resolution_clock::now() - start;
		long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

		// render grid
		ImDrawList *draw_list = ImGui::GetWindowDrawList();
		float sum = 0.0f;
		float v_sum = 0.0f;
		for (int yi = 0; yi < GRID_SIZE; yi++) {
			for (int xi = 0; xi < GRID_SIZE; xi++) {
				ImVec2 p0 = {xi * cell_size, yi * cell_size};
				ImVec2 p1 = {p0.x + cell_size, p0.y + cell_size};

				auto h = sim.water_h.at(xi, yi);
				auto u = sim.vel_u.at(xi, yi);
				auto v = sim.vel_v.at(xi, yi);
				auto uv = abs(u) + abs(v);

				auto div = sim.get_divergence(xi, yi);

				sum += h;
				v_sum += uv;

				auto color = ImColor(
					h,

					div / 5.0f,
					//uv / 10.0f,

					sim.ground_h.at(xi, yi),
					1.0f);

				draw_list->AddRectFilled(p0, p1, color);

			}
		}

		ImGui::Text("total water: %f", sum);
		ImGui::Text("total velocity: %f", v_sum);
		ImGui::Text("sim time(ms): %f", float(microseconds) * 0.001f);

		ImGui::End();

		// Rendering
		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		ImVec4 clear_color = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

	return 0;
}