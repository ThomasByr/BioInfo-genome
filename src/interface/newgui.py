import os
import logging
import re
from multiprocessing.pool import ThreadPool
import tkinter as tk
from tkinter import ttk
from CTkMessagebox import CTkMessagebox
import customtkinter as ctk
from ttkwidgets.autocomplete import AutocompleteCombobox

from PIL import Image, ImageTk

from .checkboxtreeview import CheckboxTreeview
from ..core import Tree, create_data_from_stuff

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("dark-blue")


checkboxes = {
    "-CDS-": "CDS",
    "-CENTROMERE-": "centromere",
    "-INTRON-": "intron",
    "-MOBILE_ELEMENT-": "mobile_element",
    "-NCRNA-": "ncRNA",
    "-RRNA-": "rRNA",
    "-TELOMERE-": "telomere",
    "-TRNA-": "tRNA",
    "-3UTR-": "3'UTR",
    "-5UTR-": "5'UTR",
}
depths = {
    1: "Kingdom",
    2: "Specie",
    3: "Subspecie",
    4: "Organism",
    5: "CDS Record",
}

__all__ = ["App"]

logger = logging.getLogger("newgui")


def open_file(file_path: str):
    """Opens a file using the default program for the file type"""
    if not os.path.isfile(file_path):
        return
    if os.name == "nt":
        os.startfile(file_path)
    elif os.name == "posix":
        try:
            os.system(f"wslview {file_path}")
        except:  # noqa
            try:
                os.system(f"xdg-open {file_path}")
            except:  # noqa
                raise NotImplementedError("please install xdg-open (or wslu)")  # noqa

    else:
        raise NotImplementedError(f"Unsupported OS: {os.name}")


def open_file_from_selection(selection: tuple[str, ...]):
    try:
        open_file(selection[0])
    except IndexError:
        pass


def run_internal(organism: str, selected_regions: list[str], tree: Tree) -> tuple[int, str]:
    value = tree.get_info(organism)
    return create_data_from_stuff(value.name, value.path, value.nc, selected_regions), organism


class App(ctk.CTk):
    def __init__(self, tree: Tree):
        super().__init__()
        self.__tree = tree
        self.__all_organisms = self.__tree.organisms
        self.__all_organisms_sorted = sorted(self.__all_organisms)
        self.pool = ThreadPool(processes=4)  #! 4 is max of concurrent requests on server

        # treeview Customisation (theme colors are selected)
        self.__bg_color = self._apply_appearance_mode(ctk.ThemeManager.theme["CTkFrame"]["fg_color"])
        self.__text_color = self._apply_appearance_mode(ctk.ThemeManager.theme["CTkLabel"]["text_color"])
        self.__selected_color = self._apply_appearance_mode(ctk.ThemeManager.theme["CTkButton"]["fg_color"])

        customstyle = ttk.Style()
        customstyle.theme_use("default")
        for w in (
            "Treeview",
            "TCombobox",
            "TEntry",
            "TScrollbar",
            "Treeview.Heading",
            "Treeview.Cell",
            "TLabel",
        ):
            customstyle.configure(
                w, background=self.__bg_color, foreground=self.__text_color, fieldbackground=self.__bg_color
            )
            customstyle.map(
                w,
                background=[("selected", self.__bg_color)],
                foreground=[("selected", self.__selected_color)],
            )
        self.bind("<<TreeviewSelect>>", lambda _: self.focus_set())

        # configure window
        self.title("BioInfo-Genome")
        self.geometry(f"{1100}x{580}")
        icon = ImageTk.PhotoImage(Image.open(os.path.join("assets", "favicon.png")))
        self.wm_iconbitmap()
        self.iconphoto(True, icon)

        # left frame
        self.file_tree_frame = ctk.CTkFrame(self, width=450)
        self.file_tree_frame.pack(side="left", fill="y")

        # right frame
        self.right_frame = ctk.CTkFrame(self)
        self.right_frame.pack(side="right", fill="both", expand=False)

        self.checkbox_frame = ctk.CTkFrame(self.right_frame)
        self.checkbox_frame.pack(side="top", fill="both", expand=False)

        # create the checkboxes
        # max 2 checkboxes per line
        self.checkbox_list: list[ctk.CTkCheckBox] = []
        self.checkbox_frame_list: list[ctk.CTkFrame] = []
        checkbox_per_line = 2
        for i, (_, value) in enumerate(checkboxes.items()):
            if i % checkbox_per_line == 0:
                self.checkbox_frame_list.append(ctk.CTkFrame(self.checkbox_frame))
                self.checkbox_frame_list[-1].grid(
                    row=i // checkbox_per_line, column=0, sticky="nsew", padx=5, pady=5
                )
            self.checkbox_list.append(
                ctk.CTkCheckBox(
                    self.checkbox_frame_list[-1],
                    text=value,
                    variable=tk.BooleanVar(),
                    onvalue=True,
                    offvalue=False,
                )
            )
            self.checkbox_list[-1].pack(side="left")
        self.checkbox_list[0].select()

        # preview frame
        self.preview_frame = ctk.CTkFrame(self)
        self.preview_frame.pack(side="top", fill="both", expand=True)

        # bellow, put a text box to view logs
        self.log_frame = ctk.CTkFrame(self.right_frame)
        self.log_frame.pack(side="bottom", fill="both", expand=True)
        self.log_text = ctk.CTkTextbox(self.log_frame)
        self.log_text.pack(side="bottom", fill="both", expand=True)
        self.log_text.insert("end", "Logs will appear here\n")
        self.log_text.configure(font=("TkFixedFont", 10), state="disabled")

        # create a small one-line input to search for organisms in the tree
        self.search_bar = AutocompleteCombobox(
            self.file_tree_frame, completevalues=self.__all_organisms_sorted
        )
        self.search_bar.pack(side="top", fill="x", padx=5, pady=5)
        self.search_bar.bind("<KeyRelease>", func=lambda e: self.search_bar_callback(e))
        self.search_bar.bind("<Return>", func=lambda e: self.search_bar_selection_callback(e))

        # bottom frame
        self.bottom_frame = ctk.CTkFrame(self)
        self.bottom_frame.pack(side="bottom", fill="x")

        # preview stuff
        self.preview_label = ctk.CTkLabel(self.preview_frame, text="No selection")
        self.preview_label.pack(pady=10)
        self.preview_box = ctk.CTkTextbox(self.preview_frame)
        self.preview_box.configure(state="disabled")
        self.preview_box.pack(fill="both", expand=True)

        # put a floating button on top of the preview box with 3 dashes to open a sliding menu
        self.menu_button = ctk.CTkButton(
            self.preview_box,
            text="≡",
            width=30,
            height=30,
            command=lambda: self.toggle_menu(self.tree_view.selection()),
        )
        self.menu_button.place(relx=0.01, rely=0.01)
        # color the button soft gray
        self.menu_button.configure(fg_color="#696969", hover_color="#757575")
        self.menu_button.configure(state="disabled")

        # put a floating button on top of the preview box next to the menu button to open the file
        self.open_file_button = ctk.CTkButton(
            self.preview_box,
            text="open",
            width=30,
            height=30,
            command=lambda: open_file_from_selection(self.tree_view.selection()),
        )
        self.open_file_button.place(relx=0.08, rely=0.01)
        # color the button soft gray
        self.open_file_button.configure(fg_color="#696969", hover_color="#757575")
        self.open_file_button.configure(state="disabled")

        self.__all_dynamic_buttons = [self.open_file_button, self.menu_button]

        def write_in_preview(text: str):
            self.preview_box.configure(state="normal")
            self.preview_box.delete("1.0", "end")
            self.preview_box.insert("end", text)
            self.preview_box.configure(state="disabled")

        def update_preview(selection: str):
            base_selection = os.path.basename(selection)
            if os.path.isfile(selection):
                for button in self.__all_dynamic_buttons:
                    button.configure(state="normal")
                self.preview_label.configure(text=f"Preview for {base_selection}")
                with open(selection, "r") as f:
                    write_in_preview(f.read())
            else:
                for button in self.__all_dynamic_buttons:
                    button.configure(state="disabled")
                self.preview_label.configure(text=f"Not a file: {base_selection}")
                depth = len(os.path.normpath(selection).split(os.sep)) - 1
                size = self.tree_view.item(selection)["values"][0]
                write_in_preview("")
                if 0 < depth < 5:
                    write_in_preview(
                        "\n\n\n"
                        f"{depths[depth]} {base_selection} :"
                        f" {size} {depths[depth+1]}{'s' if size > 1 else ''}",
                    )

        def on_tree_select(_: tk.Event):
            selected_item = self.tree_view.selection()[0]
            update_preview(selected_item)

        # file tree
        self.tree_scrollbar = ctk.CTkScrollbar(self.file_tree_frame)
        self.tree_scrollbar.pack(side="right", fill="y")

        self.tree_view = CheckboxTreeview(self.file_tree_frame, yscrollcommand=self.tree_scrollbar.set)
        self.tree_view.pack(expand=True, fill="both")

        self.tree_scrollbar.configure(command=self.tree_view.yview)

        self.tree_view["columns"] = ("Size",)
        self.tree_view.column("#0", width=400, minwidth=400, stretch=tk.YES)
        self.tree_view.column("Size", anchor=tk.W, width=50)
        self.tree_view.heading("#0", text="File/Folder")
        self.tree_view.heading("Size", text="Size")
        self.tree_view.bind("<<TreeviewSelect>>", on_tree_select)
        self.tree_view.bind("<Double-1>", lambda _: open_file(self.tree_view.selection()[0]))
        self.fill_tree_view("Results", "")

        # bottom frame stuff
        self.button1 = ctk.CTkButton(
            self.bottom_frame, text="run", fg_color="#556b2f", command=self.do_stuff_run
        )
        self.button1.grid(row=0, column=0, padx=5, pady=5)

        self.button2 = ctk.CTkButton(
            self.bottom_frame, text="Collapse All", command=self.tree_view.collapse_all
        )
        self.button2.grid(row=1, column=0, padx=5, pady=5)

        self.button3 = ctk.CTkButton(self.bottom_frame, text="Expand All", command=self.tree_view.expand_all)
        self.button3.grid(row=1, column=1, padx=5, pady=5)

        self.button4 = ctk.CTkButton(
            self.bottom_frame, text="Reset", fg_color="#6d3c3c", command=self.reset_tree_view
        )
        self.button4.grid(row=1, column=2, padx=5, pady=5)

        self.__all_buttons = [self.button1, self.button2, self.button3, self.button4]

        self.made_with_love_label = ctk.CTkLabel(
            self.bottom_frame, text="@ThomasByr, @m7415, @JBrandstaedt and @Bas6700 | Made with ❤️"
        )
        self.made_with_love_label.grid(row=2, column=0, columnspan=3, pady=5)

    def toggle_menu(self, selections: tuple[str, ...]):
        # show a menu with clickable elements showing lines containing:
        # "NC_016049:"
        try:
            file_path = selections[0]
        except IndexError:
            return
        if not os.path.isfile(file_path):
            return

        with open(file_path, "r") as f:
            lines = f.readlines()
        matched_lines = []
        for line in lines:
            if re.search(r"NC_\d+:", line):
                matched_lines.append(line)
        if not matched_lines:
            return

        def __callback(line: str):
            # scroll to the line
            line_number = lines.index(line) + 1
            self.preview_box.see(f"{line_number}.0")
            # unselect all
            self.preview_box.tag_remove("sel", "1.0", "end")
            # select the line
            self.preview_box.tag_add("sel", f"{line_number}.0", f"{line_number}.end")
            # focus the preview box
            self.preview_box.focus_set()

        menu = tk.Menu(self)
        menu.configure(bg=self.__bg_color, fg=self.__text_color)

        for line in matched_lines:
            menu.add_command(label=line, command=lambda line=line: __callback(line))
        menu.post(self.menu_button.winfo_rootx(), self.menu_button.winfo_rooty())

    def search_bar_callback(self, event: tk.Event):
        """
        Callback for the search bar\\
        This is called whenever the user types something in the search bar
        """
        # get the current value of the search bar
        value = event.widget.get()
        # if the value is empty, show all organisms
        if not value or value.isspace():
            self.search_bar.configure(completevalues=self.__all_organisms_sorted)
            return
        # filter out all organisms that do not match the search bar
        filtered = [x for x in self.__all_organisms if value.lower() in x.lower()]
        # update the search bar
        self.search_bar.configure(completevalues=sorted(filtered))

    def search_bar_selection_callback(self, event: tk.Event):
        """
        Callback for the search bar\\
        This is called whenever the user presses enter in the search bar
        """
        # get the current value of the search bar
        value = event.widget.get()
        # if the value is empty, show all organisms
        if not value:
            self.search_bar.configure(completevalues=self.__all_organisms_sorted)
            return
        # filter out all organisms that do not match the search bar
        filtered = [x for x in self.__all_organisms if value.lower() in x.lower()]
        __filtered = list(map(lambda x: self.__tree.get_info(x).path, filtered))
        # update the search bar
        self.search_bar.configure(completevalues=sorted(filtered))
        # if there is only one organism, select it
        if len(filtered) == 1:
            self.tree_view.focus(__filtered[0])
            self.tree_view.selection_set(__filtered[0])
            self.tree_view._check_ancestor(__filtered[0])
            self.tree_view.see(__filtered[0])
            self.tree_view.focus_set()
            self.search_bar.set("")
        # if there are multiple organisms, show a list
        elif len(filtered) > 1:
            self.search_bar.set("")
            msg = CTkMessagebox(
                master=self,
                title="Multiple organisms found",
                width=600,
                message=f"Multiple organisms found ({len(filtered)}):\n{', '.join(filtered[:5])}...",
                icon="info",
                option_1="Cancel",
                option_2="Select all",
            )
            msg.wait_window()
            if msg.get() == "Select all":
                for organism in __filtered:
                    self.tree_view.focus(organism)
                    self.tree_view.selection_set(organism)
                    self.tree_view._check_ancestor(organism)
                    self.tree_view.see(organism)
                    self.tree_view.focus_set()
            self.search_bar.set("")
        # if there are no organisms, show a warning
        else:
            self.search_bar.set("")
            msg = CTkMessagebox(
                master=self,
                title="No organisms found",
                width=600,
                message=f"No organisms found for '{value}'",
                icon="cancel",
                option_1="OK",
            )
            msg.wait_window()
            if msg.get() == "OK":
                self.search_bar.set("")

        # in all cases, reset the search bar complete values
        self.search_bar.configure(completevalues=self.__all_organisms_sorted)

    def fill_tree_view(self, folder, parent):
        """Recursively fill the tree view with the files and folders in the given folder"""
        full_path = os.path.join(parent, folder) if parent else folder
        number_of_files = len(os.listdir(full_path)) if os.path.isdir(full_path) else ""
        size = os.path.getsize(full_path) if os.path.isfile(full_path) else number_of_files

        # insert the current folder or file into the tree view
        try:
            self.tree_view.insert(parent, "end", full_path, text=folder, values=(size,))
        except tk.TclError:
            pass

        if os.path.isdir(full_path):
            for item in os.listdir(full_path):
                self.fill_tree_view(item, full_path)

    def reset_tree_view(self):
        def __callback():
            self.tree_view.delete(*self.tree_view.get_children())
            self.fill_tree_view("Results", "")
            self.__all_organisms = self.__tree.organisms
            self.__all_organisms_sorted = sorted(self.__all_organisms)
            self.search_bar.configure(completevalues=self.__all_organisms_sorted)
            self.emit_log("Reset done ✅")
            for button in self.__all_buttons:
                button.configure(require_redraw=True, state="normal")

        # show some retry/cancel warnings
        msg = CTkMessagebox(
            master=self,
            title="Warning!",
            width=600,
            message="Performing this action will reset the tree view as well as delete all the files you have downloaded.\n\nAre you sure you want to continue?",
            icon="warning",
            option_1="Cancel",
            option_2="Continue",
        )
        msg.wait_window()
        if msg.get() == "Continue":
            self.emit_log("Reset requested (please wait) ❗")
            for button in self.__all_buttons:
                button.configure(state="disabled")
            self.pool.apply_async(self.__tree.build, args=(True, True), callback=lambda _: __callback())

    def emit_log(self, text: str):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", text + "\n")
        logger.info("%s", text)
        self.log_text.configure(state="disabled")
        self.log_text.see("end")

    def do_stuff_run(self):
        def __callback(x: int, o: str):
            self.fill_tree_view("Results", "")  # todo: do not check for all tree
            self.emit_log(f"{o} {selected_regions} ({x})")

        selected_regions: list[str] = []
        for checkbox in self.checkbox_list:
            if checkbox.get():
                selected_regions.append(checkbox.cget("text"))

        selected_organisms: set[str] = set()
        print(self.tree_view.get_checked())
        for item in self.tree_view.get_checked():
            # item is a path
            # if it is a file, remove the file name
            # keep only the last folder
            if os.path.isfile(item):
                item = os.path.dirname(item)
            # get the last folder
            item = os.path.basename(item)
            selected_organisms.add(item)
        self.tree_view.uncheck_all()

        for organism in selected_organisms:
            self.pool.apply_async(
                run_internal,
                args=(organism, selected_regions, self.__tree),
                callback=lambda x: __callback(*x),
            )
