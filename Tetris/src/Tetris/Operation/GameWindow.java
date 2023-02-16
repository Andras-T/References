package Tetris.Operation;
import Tetris.Graphics.*;
import javax.swing.*;

/**
 * A játék ablaka
 */
public class GameWindow extends JFrame{

    private static final int WIDTH = 14*26+15, HEIGHT= 24*26+38;
    private final MainMenu mainMenu;
    private final Board board;
    private Game game;
    // Legmagasabb eredmények
    private HighScores scores;

    /**
     * Létrehozza az ablakot és megynitja a főmenüt
     */
    public GameWindow(){
    super("Tetris");

    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    setSize(WIDTH, HEIGHT);
    setResizable(false);
    setLocationRelativeTo(null);
    mainMenu = new MainMenu(this);
    board = new Board();
    add(mainMenu);
    setVisible(true);
    }

    /**
     *  Megnyitja a főmenüt egy másik panelről
     * @param panel lecserélendő panel
     */
    public void mainMenu(JPanel panel){
        remove(panel);
        add(mainMenu);
        setVisible(true);
        repaint();
    }

    /**
     * Elindítja a játékot
     */
    public void start(){
        remove(mainMenu);
        game = new Game(this);
        add(board);
        setVisible(true);
        board.requestFocusInWindow();
    }

    /**
     * Befejezi a játékot és megjeleníti az elért eredményt
     * @param score elért eredmény
     */
    public void gameOver(int score){
        remove(board);
        GameOver over = new GameOver(this, score);
        add(over);
        setVisible(true);
        over.requestFocusInWindow();
    }

    /**
     * Megjeleníti az eddigi 5 legmagasasbb eredményt
     */
    public void highScores(){
        remove(mainMenu);
        scores = new HighScores(this);
        add(scores);
        setVisible(true);
        scores.requestFocusInWindow();
    }

    /**
     * Megadja az aktuális játékot
     * @return  Aktuális játék
     * @throws NullPointerException ha nem indult el még játék
     */
    public Game getGame() throws NullPointerException{
        if (game == null)
            throw new NullPointerException();
        return game;
    }

    /**
     * Megadja az aktuális játékhoz tartozó panelt
     * @return aktuális játékhoz tartozó panel
     */
    public Board getBoard(){
        return board;
    }

    public HighScores getHighScores(){
        return scores;
    }

    /**
     * Bezárja az ablakot
     */
    public void close(){
        if (JOptionPane.showConfirmDialog( this,"Confirm if you want to Quit","Tetris",
                JOptionPane.YES_NO_OPTION)==JOptionPane.YES_OPTION)
        dispose();
    }

    /**
     * Megadja az ablak szélességét
     * @return ablak szélessége
     */
    public int getWidth(){
        return WIDTH;
    }

    /**
     * Megadja az ablak magasságát
     * @return ablak magassága
     */
    public int getHeight() {
        return HEIGHT;
    }
}