import os
import logging
import copy
from multiprocessing.pool import ThreadPool
import tkinter as tk
from tkinter import ttk
from CTkMessagebox import CTkMessagebox
import customtkinter as ctk

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

__all__ = ["App"]

logger = logging.getLogger("newgui")


def run_internal(organism: str, selected_regions: list[str], tree: Tree) -> tuple[int, str]:
    value = tree.get_info(organism)
    return create_data_from_stuff(value.name, value.path, value.nc, selected_regions), organism


class App(ctk.CTk):
    def __init__(self, tree: Tree):
        super().__init__()
        self.__tree = tree
        self.pool = ThreadPool(processes=4)  #! 4 is max of concurrent requests on server

        # treeview Customisation (theme colors are selected)
        bg_color = self._apply_appearance_mode(ctk.ThemeManager.theme["CTkFrame"]["fg_color"])
        text_color = self._apply_appearance_mode(ctk.ThemeManager.theme["CTkLabel"]["text_color"])
        selected_color = self._apply_appearance_mode(ctk.ThemeManager.theme["CTkButton"]["fg_color"])

        treestyle = ttk.Style()
        treestyle.theme_use("default")
        treestyle.configure(
            "Treeview",
            background=bg_color,
            foreground=text_color,
            fieldbackground=bg_color,
            borderwidth=0,
        )
        treestyle.map(
            "Treeview",
            background=[("selected", bg_color)],
            foreground=[("selected", selected_color)],
        )
        self.bind("<<TreeviewSelect>>", lambda event: self.focus_set())

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

        # bottom frame
        self.bottom_frame = ctk.CTkFrame(self)
        self.bottom_frame.pack(side="bottom", fill="x")

        # preview stuff
        self.preview_label = ctk.CTkLabel(self.preview_frame, text="No selection")
        self.preview_label.pack(pady=10)
        self.preview_box = ctk.CTkTextbox(self.preview_frame)
        self.preview_box.configure(state="disabled")
        self.preview_box.pack(fill="both", expand=True)

        def update_preview(selection):
            if os.path.isfile(selection):
                self.preview_label.configure(text=f"Preview for {os.path.basename(selection)}")
                with open(selection, "r") as f:
                    self.preview_box.configure(state="normal")
                    self.preview_box.delete("1.0", "end")
                    self.preview_box.insert("end", f.read())
                    self.preview_box.configure(state="disabled")
            else:
                self.preview_label.configure(text=f"Not a file: {os.path.basename(selection)}")
                self.preview_box.configure(state="normal")
                self.preview_box.delete("1.0", "end")
                self.preview_box.configure(state="disabled")

        def on_tree_select(_):
            selected_item = self.tree_view.selection()[0]
            update_preview(selected_item)

        # file tree
        self.tree_view = CheckboxTreeview(self.file_tree_frame)
        self.tree_view.pack(expand=True, fill="both")
        self.tree_view["columns"] = ("Size",)
        self.tree_view.column("#0", width=400, minwidth=400, stretch=tk.YES)
        self.tree_view.column("Size", anchor=tk.W, width=50)
        self.tree_view.heading("#0", text="File/Folder")
        self.tree_view.heading("Size", text="Size")
        self.tree_view.bind("<<TreeviewSelect>>", on_tree_select)
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

    def fill_tree_view(self, folder, parent):
        full_path = os.path.join(parent, folder) if parent else folder
        size = os.path.getsize(full_path) if os.path.isfile(full_path) else ""

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
            self.emit_log("Reset done ✅")
            for button in self.__all_buttons:
                button.configure(require_redraw=True, state="normal")

        # Show some retry/cancel warnings
        msg = CTkMessagebox(
            master=self,
            title="Warning!",
            width=600,
            message="Performing this action will reset the tree view as well as delete all the files you have downloaded. Are you sure you want to continue?",
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
        logger.info(text)
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
