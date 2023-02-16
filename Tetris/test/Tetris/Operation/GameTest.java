package Tetris.Operation;

import Tetris.Graphics.*;
import org.junit.Assert;
import org.junit.Test;
import java.awt.*;
import java.io.IOException;
import static org.junit.Assert.fail;

/**
 * Játékot tesztel
 */
public class GameTest {

    //GameWindow window;
    Game game;

    public GameTest(){
        game = new Game(new GameWindow());
    }

    /**
     * Tetramino elforgatása
     */
    @Test
    public void rotateTetraminoTest(){
        Tetramino original = new ShapeL();
        Tetramino rotate = original.copy();
        for (int c = 0; c < 4; c++) {
            for (int r = 0; r < 4; r++)
                Assert.assertEquals((original.getShape())[c][r], (rotate.getShape())[c][r]);
        }
        game.rotate(rotate);
        game.rotate(rotate);
        game.rotate(rotate);
        game.rotate(rotate);

        for (int c = 0; c < 4; c++) {
            for (int r = 0; r < 4; r++)
                Assert.assertEquals((original.getShape())[c][r], (rotate.getShape())[c][r]);
        }
    }

    /**
     * Alakzat elforgatása
     */
    @Test
    public void rotateShapeTest(){
        int[][] originalShape = {
                {1,1,1,0,0,0},
                {1,1,1,0,0,0},
                {1,1,1,0,0,0},
                {1,1,1,0,0,0},
                {1,1,1,0,0,0},
                {1,1,1,0,0,0}
        };
        int[][] shape = originalShape.clone();

        shape = game.rotateCW(shape);
        shape = game.rotateCW(shape);

        for (int c = 0; c < 6; c++) {
            for (int r = 0; r < 6; r++)
                Assert.assertNotEquals(originalShape[c][r], shape[c][r]);
        }
    }

    /**
     * Eredmény mentése
     */
    @Test
    public void saveTest(){
        try {
            game.save(100);
        } catch (IOException e) {
            fail("Nem sikerült megnyitni a fájlt!");
        }
    }

    /**
     * Betűtípus beolvasás
     */
    @Test
    public void gameOverFontTest(){
        GameOver gameOver = new GameOver(game.gameWindow, 0);
        try {
            gameOver.getSpecialFont();
        } catch (IOException | FontFormatException e) {
            fail("Nem sikerült betölteni a betűtípust!");
        }
    }

}
