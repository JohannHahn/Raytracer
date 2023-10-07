#include <iostream>
#include "raylib.h"


size_t window_width = 900;
size_t window_height = 600;

int main()
{

	InitWindow(window_width, window_height, "TestWindow");
	SetTargetFPS(60);
	Image img = GenImageColor(window_width, window_height, BLACK);
	Texture txt = LoadTextureFromImage(img);
	while (!WindowShouldClose()) {
		BeginDrawing();
			DrawTexture(txt, 0, 0, WHITE);
		EndDrawing();
	}
	UnloadTexture(txt);
	CloseWindow();
	return 0;
}