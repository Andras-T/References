package Tetris.Graphics;

import Tetris.Operation.GameWindow;
import org.junit.Test;

import java.awt.*;
import java.io.IOException;

import static org.junit.Assert.fail;

/**
 * Megjelenítést tesztel
 */
public class GraphicsTest {

    HighScores highScores;
    MainMenu mainMenu;

    public GraphicsTest(){
        highScores = new HighScores(new GameWindow());
        mainMenu = new MainMenu(new GameWindow());
    }

    /**
     * Eredmények beolvasása
     */
    @Test
    public void loadTest(){
        try {
            highScores.load();
        } catch (IOException e) {
            fail("Nem sikerült megnyitni a fájlt!");
        }
    }

    /**
     * Betűtípus beolvasása eredményeknél
     */
    @Test
    public void scoreFontTest() {
        try {
            highScores.getSpecialFont();
        } catch (FontFormatException | IOException e) {
            fail("Nem sikerült betölteni a betűtípust!");
        }
    }

    /**
     * Betűtípus beolvasása a főmenüben
     */
    @Test
    public void menuFontTest() {
        try {
            mainMenu.getSpecialFont();
        } catch (FontFormatException | IOException e) {
            fail("Nem sikerült betölteni a betűtípust!");
        }
    }
}